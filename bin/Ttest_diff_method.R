
#######################################################################
# This file is part of Dekupl. 
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Jérôme Audoux (jerome.audoux@inserm.fr), Copyright (C) 2017  
#######################################################################

Ttest                     = snakemake@input$binary
no_GENCODE                = snakemake@input$counts
normalization_factor_path = snakemake@input$sample_conditions
pvalue_threshold          = snakemake@params$pvalue_threshold
log2fc_threshold          = snakemake@params$log2fc_threshold
conditionA                = snakemake@params$conditionA
conditionB                = snakemake@params$conditionB

output_diff_counts  = snakemake@output$diff_counts
output_pvalue_all   = snakemake@output$pvalue_all
output_tmp          = snakemake@output$tmp_dir
output_log          = snakemake@log[[1]]

# Function for logging to the output
logging <- function(str) {
  sink(file=paste(output_log), append=TRUE, split=TRUE)
  print(paste(Sys.time(),str))
  sink()
}

dir.create(output_tmp, showWarnings = FALSE)

logging(paste("Start Ttest_diff_methods"))

system(paste(Ttest ,"-p", pvalue_threshold, "-f", log2fc_threshold, "-r", output_pvalue_all, no_GENCODE, normalization_factor_path, conditionA, conditionB, "| gzip -c > ",output_diff_counts))

logging(paste("End Ttest_diff_methods"))
