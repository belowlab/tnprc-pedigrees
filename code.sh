#variant filtering for full coverage data:
input_vcf=starting_vcf_file.vcf
gatk VariantFiltration \
   --variant ${input_vcf} \
   --output tulane_allchr_filtered.vcf \
   --filter-name "hard_filter" \
   --filter-expression "AC == 0 || QD < 10.0 || FS > 2.0 || ( MQ < 59.0 && MQ > 61.00 ) || SOR > 1.5"

#variant filtering for downsampled data
bcftools filter -i 'FORMAT/DP>5 && GQ>12' -S . -o downsampled.5.0.extfilt.vcf starting_vcf_file.vcf

#running PRIMUS:
run_PRIMUS.pl --file downsampled.5.0.extfilt --genome --keep_inter_files --internal_ref --output_dir filt0217_PRIMUS

#splitting by chromosome:
for c in {1..20}; do plink --bfile downsampled.5.0.qualfilt.lowmiss --chr ${c} --recode vcf --out downsampled.5.0.filtered.chr${c}; done

#phasing:
for i in {1..20}; do echo "java -jar beagle.25Nov19.28d.jar gt=downsampled.5.0.filtered.chr${i}.vcf nthreads=3 out=downsampled.5.0.chr${i}.phased"; done | parallel

#converting phased to PLINK format:
vcftools --gzvcf downsampled.5.0.chr4.phased.vcf.gz --plink --out downsampled.5.0.chr4.phased.plink

#adding position in CM:
awk '{ OFS="\t"; $3=($4/1000000)*0.433; print; }' downsampled.5.0.chr4.phased.plink.map > downsampled.5.0.chr4.phased.plink.cm.map

#running GERMLINE:
het=1
hom=2
for c in {1..20}
do
	echo "1" > germline_chr${c}.run
	echo "downsampled.5.0.ersa.chr${c}.phased.plink.cm.map" >> germline_chr${c}.run
	echo "downsampled.5.0.ersa.chr${c}.phased.plink.ped" >> germline_chr${c}.run
	echo "downsampled.5.0.ersa.chr${c}.phased.germline" >> germline_chr${c}.run
done
for c in {1..20}
do
	echo "germline -min_m 2.5 -err_het ${het} -err_hom ${hom} < germline_chr${c}.run"
done | parallel

#running ERSA
ersa --segment_files=*.match --number_of_chromosomes=20 --rec_per_meioses=13.6239623 --confidence_level=0.999 --output_file=ersa_results.txt --model_output_file=ersa_models.txt

#varying recombinations per meiosis in ERSA:
for r in $(seq 6.0 0.2 13.4)
do
        ersa --segment_files=*.match --number_of_chromosomes=20 --confidence_level=0.999 --control_files=nontulane_0419.match --mask_common_shared_regions=true --mask_region_threshold=6 --rec_per_meioses=${r} --output_file=ersa_results_varied_rpm_${r}.txt --model_output_file=ersa_models_varied_rpm_${r}.txt
done

#running PADRE
run_PRIMUS.pl --project_summary downsampled.5.0_PRIMUS/Summary_downsampled.5.0.genome.txt --ersa_model_output ersa_models.txt --ersa_results ersa_results.txt --degree_rel_cutoff 1
