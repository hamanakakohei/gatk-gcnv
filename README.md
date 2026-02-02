# gatk-gcnv

# 注意点（色々なページから集めた）
1. https://github.com/broadinstitute/gatk/tree/master/scripts/cnv_wdl/germline
CNVGermlineCohortWorkflow.intervals -- Picard or GATK-style interval list. For WGS, this should typically only include the chromosomes of interest.
CNVGermlineCohortWorkflow.num_intervals_per_scatter -- Number of intervals (i.e., targets or bins) in each scatter for GermlineCNVCaller. If total number of intervals is not divisible by the value provided, the last scatter will contain the remainder.
CNVGermlineCohortWorkflow.maximum_number_events_per_sample -- Maximum number of events threshold for doing sample QC (recommended for WES is ~100)
CNVGermlineCohortWorkflow.do_explicit_gc_correction -- (optional) If true, perform explicit GC-bias correction when creating PoN and in subsequent denoising of case samples. If false, rely on PCA-based denoising to correct for GC bias.
CNVGermlineCohortWorkflow.PreprocessIntervals.bin_length -- Size of bins (in bp) for coverage collection. This must be the same value used for all case samples.
CNVGermlineCohortWorkflow.PreprocessIntervals.padding -- Amount of padding (in bp) to add to both sides of targets for WES coverage collection. This must be the same value used for all case samples.

2. https://gatk.broadinstitute.org/hc/en-us/articles/360041850711-GermlineCNVCaller
The quality of coverage model parametrization and the sensitivity/precision of germline CNV calling are sensitive to the choice of model hyperparameters, including
the prior probability of alternative copy-number states (set using the p-alt argument)
prevalence of active (i.e. CNV-rich) intervals (set via the p-active argument)
the coherence length of CNV events and active/silent domains across intervals (set using the cnv-coherence-length and class-coherence-length arguments, respectively)
the coherence length of CNV events and active/silent domains across intervals (set using the cnv-coherence-length and class-coherence-length arguments, respectively)
the typical scale of interval-specific unexplained variance (set using the interval-psi-scale).
the typical scale of sample-specific unexplained variance (set using sample-psi-scale arguments).
It is crucial to note that these hyperparameters are not universal and must be tuned for each sequencing protocol and properly set at runtime.

 For WES and WGS, we recommend no less than 10000 consecutive intervals spanning at least 10 - 50 mb.

Memory requirements for the python subprocess ("gcnvkernel"):
The computation done by this tool, for the most part, is performed outside of JVM and via a spawned python subprocess. The Java heap memory is only used for loading sample counts and preparing raw data for the python subprocess.
The user must ensure that the machine has enough free physical memory for spawning and executing the python subprocess.
Generally speaking, the resource requirements of this tool scale linearly with each of
the number of samples
the number of modeled intervals
the highest copy number state
the number of bias factors
the number of knobs on the GC curve.
For example, the python subprocess requires approximately 16GB of physical memory for modeling
10000 intervals for 100 samples
with 16 maximum bias factors
maximum copy-number state of 10
explicit GC bias modeling.

wgs
paddingを0にするのを忘れずに
bin lengthはデフォの1000
do_explicit_gc_correctionはデフォでtrue
mappability_track_bedはhoffmanラボのデータを
segmental_duplication_track_bedこっちはどこから？
blacklist regionsでセントロメアを指定？
cnv_coherence_length 1000
class_coherence_length 1000
enable_bias_factors false
interval-psi-scale 1e-6
log-mean-bias-standard-deviation 0.01
sample-psi-scale 1e-6
--> QS filtering?

__WGS parameters that increase the sensitivity of calling from Mark Walker
    --class-coherence-length 1000.0 \
    --cnv-coherence-length 1000.0 \
    --enable-bias-factors false \
    --interval-psi-scale 1.0E-6 \
    --log-mean-bias-standard-deviation 0.01 \
    --sample-psi-scale 1.0E-6 \

Comments on select sensitivity parameters

Decreasing --class-coherence-length from its default of 10,000bp to 1000bp decreases the expected length of contiguous segments. Factor for bin size when tuning.
Decreasing --cnv-coherence-length from its default 10,000bp to 1000bp decreases the expected length of CNV events. Factor for bin size when tuning.
Turning off --enable-bias-factors from the default true state to false turns off active discovery of learnable bias factors. This should always be on for targeted exome data and in general can be turned off for WGS data.
Decreasing --interval-psi-scale from its default of 0.001 to 1.0E-6 reduces the scale the tool considers normal in per-interval noise.
Decreasing --log-mean-bias-standard-deviation from its default of 0.1 to 0.01 reduces what is considered normal noise in bias factors.
Decreasing --sample-psi-scale from its default of 0.0001 to 1.0E-6 reduces the scale that is considered normal in sample-to-sample variance.
Additional parameters to consider include --depth-correction-tau, --p-active and --p-alt.

--depth-correction-tau has a default of 10000.0 (10K) and defines the precision of read-depth concordance with the global depth value.
--p-active has a default of 1e-2 (0.01) and defines the prior probability of common CNV states.
p-alt has a default of 1e-6 (0.000001) and defines the expected probability of CNV events (in rare CNV states).


We normalized the coverage profile by the median
sample coverage, performed PCA dimensionality reduction, and clustered samples based on the first
three principal components.


wes
paddingそのままの250でよい
bin lengthは0
-Lでキャプチャー領域を指定
リードカウントによるインターバルのフィルターはこのままでよいのだろうか、、、
enable_bias_factors true

gatk annotateintervals時にいろいろやる：
The tool requires the -R reference and the -L intervals. The tool calculates GC-content for the intervals using the reference.Although optional for the tool, we recommend annotating mappability by providing a --mappability-track regions file in either .bed or .bed.gz format. Be sure to merge any overlapping intervals beforehand. The tutorial omits use of this resource.
GATK recommends use of the the single-read mappability track, as the multi-read track requires much longer times to process. For example, the Hoffman lab at the University of Toronto provides human and mouse mappability BED files for various kmer lengths at https://bismap.hoffmanlab.org. The accompanying publication is titled Umap and Bismap: quantifying genome and methylome mappability.
Optionally and additionally, annotate segmental duplication content by providing a --segmental-duplication-track regions file in either .bed or .bed.gz format.
Exclude undesirable intervals with the -XL parameter, e.g. intervals corresponding to centromeric regions.

Scatternoサイズについて
4. Call copy number variants with GermlineCNVCaller
GermlineCNVCaller learns a denoising model per scattered shard while consistently calling CNVs across the shards. The tool models systematic biases and CNVs simultaneously, which allows for sensitive detection of both rare and common CNVs. For a description of innovations, see Blog #23439 As the tool index states under Important Remarks (v4.1.0.0), the tool should see data from a large enough genomic region so as to be exposed to diverse genomic features. The current recommendation is to provide at least ~10–50Mbp genomic coverage per scatter. This applies to exomes or WGS. This allows reliable inference of bias factors including GC bias. The limitation of analyzing larger regions is available memory. As an analysis covers more data, memory requirements increase.

The default --max-copy-number is capped at 5. This means the tool reports any events with more copies as CN5.
