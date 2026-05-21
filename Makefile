# Directories
DATA_DIR  := outputs/data
FIGS_DIR  := outputs/figs
JULIA_DIR := Julia

# Targets
STABILITY_DATA_FILES     := $(DATA_DIR)/fig_3_data_julia.csv $(DATA_DIR)/fig_3_d_vals.csv $(DATA_DIR)/fig_3_sigma_vals.csv
STABILITY_PLOT_FILE      := $(FIGS_DIR)/fig_3_repro_julia.pdf
DOSE_RESPONSE_DATA_FILES := $(DATA_DIR)/dose-response_neg.csv $(DATA_DIR)/dose-response_null.csv $(DATA_DIR)/dose-response_pos.csv
DOSE_RESPONSE_PLOT_FILE  := $(FIGS_DIR)/background_dose_response_julia.pdf

.PHONY: all compute_stability plot_stability compute_dose_response plot_dose_response clean_stability clean_dose_response clean

all: plot_stability plot_dose_response

compute_stability: $(STABILITY_DATA_FILES)

$(STABILITY_DATA_FILES): $(JULIA_DIR)/parallel_stab.jl | $(DATA_DIR)
	cd $(JULIA_DIR) && julia parallel_stab.jl

plot_stability: $(STABILITY_PLOT_FILE)

$(STABILITY_PLOT_FILE): $(STABILITY_DATA_FILES) $(JULIA_DIR)/plot_stab.jl | $(FIGS_DIR)
	cd $(JULIA_DIR) && julia plot_stab.jl

compute_dose_response: $(DOSE_RESPONSE_DATA_FILES)

$(DOSE_RESPONSE_DATA_FILES): $(JULIA_DIR)/dose_response.jl | $(DATA_DIR)
	cd $(JULIA_DIR) && julia dose_response.jl

plot_dose_response: $(DOSE_RESPONSE_PLOT_FILE)

$(DOSE_RESPONSE_PLOT_FILE): $(STABILITY_DATA_FILES) $(DOSE_RESPONSE_DATA_FILES) $(JULIA_DIR)/plot_dose_response.jl | $(FIGS_DIR)
	cd $(JULIA_DIR) && julia plot_dose_response.jl

# Create output directories if they don't exist
$(DATA_DIR):
	mkdir -p $(DATA_DIR)

$(FIGS_DIR):
	mkdir -p $(FIGS_DIR)

clean_stability:
	rm -f $(STABILITY_DATA_FILES) $(STABILITY_PLOT_FILE)

clean_dose_response:
	rm -f $(DOSE_RESPONSE_DATA_FILES) $(DOSE_RESPONSE_PLOT_FILE)

clean: clean_stability clean_dose_response