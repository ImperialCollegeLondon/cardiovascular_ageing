# Requires GNU Make Version >= 4.3 (due to grouped target operator &:)
RCMD:= /opt/R/4.0.1/bin/Rscript
MAKE2GRAPHCMD:= ~/Downloads/makefile2graph/make2graph

RAWFILE:= data/raw_df.csv.gz

SRDUFILE:=data/Derived/self_reported_drugs_used.csv
DFFILE:=data/Derived/dataframes_drugs_self_reported.csv.gz

FITFILES:= data/Derived/fits.csv data/Derived/fits.glm.Rds

.PHONY: all clean

all: results/report.html results/report.tex

$(RAWFILE): code/get_data.R data/coding4.tsv data/annotated_primary_care_prescriptions_cmri_34k_cohort.csv.gz
	$(RCMD) code/get_data.R

$(DFFILE) $(SRDUFILE): $(RAWFILE) code/prepare_data.R code/functions.R
	$(RCMD) code/prepare_data.R

$(FITFILES): $(DFFILE) code/analyses.R code/functions.R
	$(RCMD) code/analyses.R

results/report.html:  $(FITFILES) $(DFFILE) $(SRDUFILE) code/report.Rmd
	$(RCMD) -e "rmarkdown::render('code/report.Rmd',output_file='../results/report.html')"

results/report.tex:  $(FITFILES) $(DFFILE) $(SRDUFILE) code/report.Rmd
	$(RCMD) -e "rmarkdown::render('code/report.Rmd',rmarkdown::latex_document(fig_crop =FALSE), output_file='../results/report.tex')"

clean:
	rm -f $(FITFILES)
	rm -f $(DFFILE)
	rm -f $(RAWFILE)
	rm -f $(SRDUFILE)
	rm -f data/*_unselected*.csv
	rm -f data/eid_andcardiacAGE_May2022_primaryphenotype.csv
	rm -f data/eid_alcohol.csv data/smoke_exercise.csv data/500k_telomere.csv data/Riskfactors_README.txt
	#rm -f dependencies.png

dependencies.png:
	make -Bnd | $(MAKE2GRAPHCMD) | dot -Tpng -o dependencies.png
