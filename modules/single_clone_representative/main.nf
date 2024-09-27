

process single_clone_representative {

	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_clone_rep-passed.tsv.*$/) "clones/${name}_clone_rep-passed.tsv"}
	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*txt$/) "clones/${name}_clone_report.tsv"}
	
	input:
		path(airrFile)

	output:
		path("*_clone_rep-passed.tsv*"), emit: output
		path("*txt")
		val(true), emit: ready

	script:
		name = params.sample_name
		outname = airrFile.toString() - '.tsv' +"_clone_rep-passed"
		outfile = outname + ".tsv"

		"""
		#!/usr/bin/env Rscript

		## functions
		# find the different position between sequences

		src <- 
		"#include <Rcpp.h>
		using namespace Rcpp;
		#include <iostream>
		#include <vector>
		#include <string>
		#include <algorithm>
		#include <unordered_set>

		// [[Rcpp::export]]

		int allele_diff(std::vector<std::string> germs) {
			std::vector<std::vector<char>> germs_m;
			for (const std::string& germ : germs) {
				germs_m.push_back(std::vector<char>(germ.begin(), germ.end()));
			}

			int max_length = 0;
			for (const auto& germ : germs_m) {
				max_length = std::max(max_length, static_cast<int>(germ.size()));
			}

			for (auto& germ : germs_m) {
				germ.resize(max_length, '.'); // Pad with '.' to make all germs equal length
			}

			auto setdiff_mat = [](const std::vector<char>& x) -> int {
				std::unordered_set<char> unique_chars(x.begin(), x.end());
				std::unordered_set<char> filter_chars = { '.', 'N', '-' };
				int diff_count = 0;
				for (const char& c : unique_chars) {
					if (filter_chars.find(c) == filter_chars.end()) {
						diff_count++;
					}
				}
				return diff_count;
			};

			std::vector<int> idx;
			for (int i = 0; i < max_length; i++) {
				std::vector<char> column_chars;
				for (const auto& germ : germs_m) {
					column_chars.push_back(germ[i]);
				}
				int diff_count = setdiff_mat(column_chars);
				if (diff_count > 1) {
					idx.push_back(i);
				}
			}

			return idx.size();
		}"

		## libraries
		require(dplyr)
		library(Rcpp)
		library(ggplot2)
		sourceCpp(code = src)

		data <- readr::read_tsv("${airrFile}")

		nreads <- nrow(data)

		# calculating mutation between IMGT sequence and the germline sequence, selecting a single sequence to each clone with the fewest mutations
		data[["mut"]] <- sapply(1:nrow(data),function(j){
			x <- c(data[['sequence_alignment']][j], data[['germline_alignment_d_mask']][j])
			allele_diff(x)
		})
		# filter to the fewest mutations
		data <- data %>% dplyr::group_by(clone_id) %>% 
					dplyr::mutate(clone_size = n())

		data <- data %>% dplyr::group_by(clone_id) %>% dplyr::slice(which.min(mut))
		cat(paste0('Note: ', nrow(data),' sequences after selecting a single representative'))
		readr::write_tsv(data, file = "${outfile}")

		lines <- c(
			paste("START>", "Selecting clonal representative"),
			paste("PASS>", nrow(data)),
			paste("FAIL>", nreads-nrow(data)),
			paste("END>", "Selecting clonal representative"),
			"",
			""
		)

		file_path <- paste("${outname}","output.txt", sep="-")

		cat(lines, sep = "\n", file = file_path, append = TRUE)

		cat(lines, sep = "\n")
		"""
}
