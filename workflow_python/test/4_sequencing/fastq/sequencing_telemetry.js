[
    {
        "aggregation": "segment",
        "analysis_id": "7f6ead3c-1603-4573-84a0-b1bba7434f92",
        "basecall_1d": {
            "exit_status_dist": {
                "fail:qscore_filter": 10
            },
            "qscore_dist_temp": [
                {
                    "count": 3,
                    "mean_qscore": 5.0
                },
                {
                    "count": 4,
                    "mean_qscore": 5.5
                },
                {
                    "count": 3,
                    "mean_qscore": 6.0
                }
            ],
            "qscore_sum_temp": {
                "count": 10,
                "mean": 5.7552285194397,
                "sum": 57.5522842407227
            },
            "read_len_events_sum_temp": 2130,
            "seq_len_bases_dist_temp": [
                {
                    "count": 10,
                    "length": 0.0
                }
            ],
            "seq_len_bases_sum_temp": 10,
            "seq_len_events_dist_temp": [
                {
                    "count": 10,
                    "length": 0.0
                }
            ],
            "speed_bases_per_second_dist_temp": [
                {
                    "count": 10,
                    "speed": 3.0
                }
            ],
            "strand_median_pa": {
                "count": 10,
                "mean": 95.6615905761719,
                "sum": 956.615905761719
            },
            "strand_sd_pa": {
                "count": 10,
                "mean": 8.95313549041748,
                "sum": 89.5313568115234
            }
        },
        "channel_count": 1,
        "context_tags": {
            "experiment_kit": "genomic_dna",
            "experiment_type": "customer_qc",
            "filename": "pc_kl_17324_20170702_fnfaf03884_mn20321_sequencing_run_mt_2_92553",
            "local_bc_comp_model": "complement_r9.4_250bps_5mer",
            "local_bc_temp_model": "template_r9.4_450bps_5mer",
            "sample_frequency": "4000",
            "user_filename_input": "mt-2"
        },
        "latest_run_time": 6340.27490234375,
        "levels_sums": {
            "count": 10,
            "mean": 203.743789672852,
            "open_pore_level_sum": 2037.43786621094
        },
        "opts": {
            "adapter_pt_range_scale": "5.200000",
            "additional_context_bases": "2",
            "align_ref": "",
            "allow_inferior_barcodes": "0",
            "arrangements_files": "barcode_arrs_pcr12.cfg barcode_arrs_pcr96.cfg barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg barcode_arrs_nb96.cfg barcode_arrs_rbk.cfg barcode_arrs_lwb.cfg barcode_arrs_rlb.cfg barcode_arrs_rab.cfg barcode_arrs_rbk4.cfg barcode_arrs_rbk096.cfg barcode_arrs_vmk.cfg barcode_arrs_vmk2.cfg barcode_arrs_16s.cfg",
            "as_cpu_threads_per_scaler": "2",
            "as_gpu_runners_per_device": "2",
            "as_model_file": "adapter_scaling_dna_r9.4.1_min.jsn",
            "as_num_scalers": "1",
            "as_reads_per_runner": "32",
            "bam_methylation_threshold": "5.000000",
            "bam_out": "0",
            "barcode_kits": "",
            "bed_file": "",
            "builtin_scripts": "1",
            "calib_detect": "0",
            "calib_max_sequence_length": "3800",
            "calib_min_coverage": "0.600000",
            "calib_min_sequence_length": "3000",
            "calib_reference": "lambda_3.6kb.fasta",
            "chunk_size": "2000",
            "chunks_per_caller": "10000",
            "chunks_per_runner": "512",
            "client_id": "-1",
            "compress_fastq": "0",
            "cpu_threads_per_caller": "8",
            "detect_mid_strand_barcodes": "0",
            "device": "",
            "disable_pings": "0",
            "dmean_threshold": "1.000000",
            "dmean_win_size": "2",
            "end_gap1": "40",
            "end_gap2": "40",
            "extend_gap1": "40",
            "extend_gap2": "160",
            "fast5_out": "0",
            "flowcell": "",
            "front_window_size": "150",
            "gpu_runners_per_device": "4",
            "high_priority_threshold": "10",
            "input_file_list": "",
            "jump_threshold": "1.000000",
            "kernel_path": "",
            "kit": "",
            "lamp_arrangements_files": "barcode_arrs_lamp8.cfg barcode_arrs_lamp96.cfg",
            "lamp_kit": "",
            "log_speed_frequency": "0",
            "max_queued_reads": "2000",
            "max_search_len": "1000",
            "medium_priority_threshold": "4",
            "min_length_lamp_context": "40.000000",
            "min_length_lamp_target": "80.000000",
            "min_qscore": "7.000000",
            "min_score": "60.000000",
            "min_score_lamp": "80.000000",
            "min_score_lamp_mask": "50.000000",
            "min_score_lamp_target": "75.000000",
            "min_score_mask": "-1.000000",
            "min_score_mid_barcodes": "40.000000",
            "min_score_rear_override": "60.000000",
            "model_file": "template_r9.4.1_450bps_hac.jsn",
            "nested_output_folder": "0",
            "num_alignment_threads": "4",
            "num_barcode_threads": "4",
            "num_barcoding_buffers": "96",
            "num_callers": "1",
            "num_extra_bases_trim": "0",
            "open_gap1": "40",
            "open_gap2": "160",
            "overlap": "50",
            "override_scaling": "0",
            "ping_segment_duration": "60",
            "ping_url": "https://ping.oxfordnanoportal.com/basecall",
            "post_out": "0",
            "print_workflows": "0",
            "progress_stats_frequency": "-1.000000",
            "pt_median_offset": "2.500000",
            "pt_minimum_read_start_index": "30",
            "pt_required_adapter_drop": "30.000000",
            "pt_scaling": "0",
            "qscore_filtering": "0",
            "qscore_offset": "0.421000",
            "qscore_scale": "0.879000",
            "quiet": "0",
            "read_batch_size": "4000",
            "read_id_list": "",
            "rear_window_size": "150",
            "records_per_fastq": "4000",
            "recursive": "1",
            "require_barcodes_both_ends": "0",
            "resume": "0",
            "reverse_sequence": "0",
            "scaling_mad": "1.000000",
            "scaling_med": "0.000000",
            "score_matrix_filename": "5x5_mismatch_matrix.txt",
            "start_gap1": "40",
            "start_gap2": "40",
            "stay_penalty": "1.000000",
            "temp_bias": "1.000000",
            "temp_weight": "1.000000",
            "trace_categories_logs": "",
            "trace_domains_config": "",
            "trace_domains_log": "",
            "trim_barcodes": "0",
            "trim_min_events": "3",
            "trim_strategy": "dna",
            "trim_threshold": "2.500000",
            "u_substitution": "0",
            "verbose_logs": "0"
        },
        "read_count": 10,
        "reads_per_channel_dist": [
            {
                "channel": 100,
                "count": 10
            }
        ],
        "run_id": "c2d19c211888bc09d8e077df271f325c911c1010",
        "segment_duration": 60,
        "segment_number": 2,
        "segment_type": "guppy-acquisition",
        "software": {
            "analysis": "1d_basecalling",
            "name": "guppy-basecalling",
            "version": "4.2.2+effbaf8"
        },
        "tracking_id": {
            "asic_id": "83717495",
            "asic_id_eeprom": "2085578",
            "asic_temp": "31.5119057",
            "auto_update": "1",
            "auto_update_source": "https://mirror.oxfordnanoportal.com/software/MinKNOW/",
            "bream_core_version": "1.5.12.1",
            "bream_is_standard": "1",
            "bream_nc_version": "1.5.12.1",
            "bream_ont_version": "1.5.12.1",
            "bream_prod_version": "1.5.12.1",
            "bream_rnd_version": "0.0.0.0",
            "device_id": "MN20321",
            "exp_script_name": "C:\\Program Files\\OxfordNanopore\\MinKNOW\\python\\recipes\\nc\\NC_48Hr_Sequencing_Run_FLO-MIN106_SQK-LSK108_plus_Basecaller.py",
            "exp_script_purpose": "sequencing_run",
            "exp_start_time": "2017-07-02T17:19:26Z",
            "flow_cell_id": "FAF03884",
            "heatsink_temp": "34.4726563",
            "hostname": "PC-KL-17324",
            "installation_type": "map",
            "local_firmware_file": "0",
            "msg_id": "68ddb295-7a9b-4369-8f33-ea04169f4a8d",
            "operating_system": "Windows 6.2",
            "protocol_run_id": "a413384e-41eb-4264-bacc-77d95cd96685",
            "protocols_version": "1.5.12",
            "run_id": "c2d19c211888bc09d8e077df271f325c911c1010",
            "sample_id": "mt-2",
            "time_stamp": "2020-10-28T09:27:01Z",
            "usb_config": "1.1.1_ONT#MinION_fpga_1.1.0#ctrl#Auto",
            "version": "1.5.12"
        }
    },
    {
        "aggregation": "cumulative",
        "analysis_id": "7f6ead3c-1603-4573-84a0-b1bba7434f92",
        "basecall_1d": {
            "exit_status_dist": {
                "fail:qscore_filter": 10
            },
            "qscore_dist_temp": [
                {
                    "count": 3,
                    "mean_qscore": 5.0
                },
                {
                    "count": 4,
                    "mean_qscore": 5.5
                },
                {
                    "count": 3,
                    "mean_qscore": 6.0
                }
            ],
            "qscore_sum_temp": {
                "count": 10,
                "mean": 5.7552285194397,
                "sum": 57.5522842407227
            },
            "read_len_events_sum_temp": 2130,
            "seq_len_bases_dist_temp": [
                {
                    "count": 10,
                    "length": 0.0
                }
            ],
            "seq_len_bases_sum_temp": 10,
            "seq_len_events_dist_temp": [
                {
                    "count": 10,
                    "length": 0.0
                }
            ],
            "speed_bases_per_second_dist_temp": [
                {
                    "count": 10,
                    "speed": 3.0
                }
            ],
            "strand_median_pa": {
                "count": 10,
                "mean": 95.6615905761719,
                "sum": 956.615905761719
            },
            "strand_sd_pa": {
                "count": 10,
                "mean": 8.95313549041748,
                "sum": 89.5313568115234
            }
        },
        "channel_count": 1,
        "context_tags": {
            "experiment_kit": "genomic_dna",
            "experiment_type": "customer_qc",
            "filename": "pc_kl_17324_20170702_fnfaf03884_mn20321_sequencing_run_mt_2_92553",
            "local_bc_comp_model": "complement_r9.4_250bps_5mer",
            "local_bc_temp_model": "template_r9.4_450bps_5mer",
            "sample_frequency": "4000",
            "user_filename_input": "mt-2"
        },
        "latest_run_time": 6340.27490234375,
        "levels_sums": {
            "count": 10,
            "mean": 203.743789672852,
            "open_pore_level_sum": 2037.43786621094
        },
        "opts": {
            "adapter_pt_range_scale": "5.200000",
            "additional_context_bases": "2",
            "align_ref": "",
            "allow_inferior_barcodes": "0",
            "arrangements_files": "barcode_arrs_pcr12.cfg barcode_arrs_pcr96.cfg barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg barcode_arrs_nb96.cfg barcode_arrs_rbk.cfg barcode_arrs_lwb.cfg barcode_arrs_rlb.cfg barcode_arrs_rab.cfg barcode_arrs_rbk4.cfg barcode_arrs_rbk096.cfg barcode_arrs_vmk.cfg barcode_arrs_vmk2.cfg barcode_arrs_16s.cfg",
            "as_cpu_threads_per_scaler": "2",
            "as_gpu_runners_per_device": "2",
            "as_model_file": "adapter_scaling_dna_r9.4.1_min.jsn",
            "as_num_scalers": "1",
            "as_reads_per_runner": "32",
            "bam_methylation_threshold": "5.000000",
            "bam_out": "0",
            "barcode_kits": "",
            "bed_file": "",
            "builtin_scripts": "1",
            "calib_detect": "0",
            "calib_max_sequence_length": "3800",
            "calib_min_coverage": "0.600000",
            "calib_min_sequence_length": "3000",
            "calib_reference": "lambda_3.6kb.fasta",
            "chunk_size": "2000",
            "chunks_per_caller": "10000",
            "chunks_per_runner": "512",
            "client_id": "-1",
            "compress_fastq": "0",
            "cpu_threads_per_caller": "8",
            "detect_mid_strand_barcodes": "0",
            "device": "",
            "disable_pings": "0",
            "dmean_threshold": "1.000000",
            "dmean_win_size": "2",
            "end_gap1": "40",
            "end_gap2": "40",
            "extend_gap1": "40",
            "extend_gap2": "160",
            "fast5_out": "0",
            "flowcell": "",
            "front_window_size": "150",
            "gpu_runners_per_device": "4",
            "high_priority_threshold": "10",
            "input_file_list": "",
            "jump_threshold": "1.000000",
            "kernel_path": "",
            "kit": "",
            "lamp_arrangements_files": "barcode_arrs_lamp8.cfg barcode_arrs_lamp96.cfg",
            "lamp_kit": "",
            "log_speed_frequency": "0",
            "max_queued_reads": "2000",
            "max_search_len": "1000",
            "medium_priority_threshold": "4",
            "min_length_lamp_context": "40.000000",
            "min_length_lamp_target": "80.000000",
            "min_qscore": "7.000000",
            "min_score": "60.000000",
            "min_score_lamp": "80.000000",
            "min_score_lamp_mask": "50.000000",
            "min_score_lamp_target": "75.000000",
            "min_score_mask": "-1.000000",
            "min_score_mid_barcodes": "40.000000",
            "min_score_rear_override": "60.000000",
            "model_file": "template_r9.4.1_450bps_hac.jsn",
            "nested_output_folder": "0",
            "num_alignment_threads": "4",
            "num_barcode_threads": "4",
            "num_barcoding_buffers": "96",
            "num_callers": "1",
            "num_extra_bases_trim": "0",
            "open_gap1": "40",
            "open_gap2": "160",
            "overlap": "50",
            "override_scaling": "0",
            "ping_segment_duration": "60",
            "ping_url": "https://ping.oxfordnanoportal.com/basecall",
            "post_out": "0",
            "print_workflows": "0",
            "progress_stats_frequency": "-1.000000",
            "pt_median_offset": "2.500000",
            "pt_minimum_read_start_index": "30",
            "pt_required_adapter_drop": "30.000000",
            "pt_scaling": "0",
            "qscore_filtering": "0",
            "qscore_offset": "0.421000",
            "qscore_scale": "0.879000",
            "quiet": "0",
            "read_batch_size": "4000",
            "read_id_list": "",
            "rear_window_size": "150",
            "records_per_fastq": "4000",
            "recursive": "1",
            "require_barcodes_both_ends": "0",
            "resume": "0",
            "reverse_sequence": "0",
            "scaling_mad": "1.000000",
            "scaling_med": "0.000000",
            "score_matrix_filename": "5x5_mismatch_matrix.txt",
            "start_gap1": "40",
            "start_gap2": "40",
            "stay_penalty": "1.000000",
            "temp_bias": "1.000000",
            "temp_weight": "1.000000",
            "trace_categories_logs": "",
            "trace_domains_config": "",
            "trace_domains_log": "",
            "trim_barcodes": "0",
            "trim_min_events": "3",
            "trim_strategy": "dna",
            "trim_threshold": "2.500000",
            "u_substitution": "0",
            "verbose_logs": "0"
        },
        "read_count": 10,
        "reads_per_channel_dist": [
            {
                "channel": 100,
                "count": 10
            }
        ],
        "run_id": "c2d19c211888bc09d8e077df271f325c911c1010",
        "segment_duration": 120,
        "segment_number": 1,
        "segment_type": "guppy-acquisition",
        "software": {
            "analysis": "1d_basecalling",
            "name": "guppy-basecalling",
            "version": "4.2.2+effbaf8"
        },
        "tracking_id": {
            "asic_id": "83717495",
            "asic_id_eeprom": "2085578",
            "asic_temp": "31.5119057",
            "auto_update": "1",
            "auto_update_source": "https://mirror.oxfordnanoportal.com/software/MinKNOW/",
            "bream_core_version": "1.5.12.1",
            "bream_is_standard": "1",
            "bream_nc_version": "1.5.12.1",
            "bream_ont_version": "1.5.12.1",
            "bream_prod_version": "1.5.12.1",
            "bream_rnd_version": "0.0.0.0",
            "device_id": "MN20321",
            "exp_script_name": "C:\\Program Files\\OxfordNanopore\\MinKNOW\\python\\recipes\\nc\\NC_48Hr_Sequencing_Run_FLO-MIN106_SQK-LSK108_plus_Basecaller.py",
            "exp_script_purpose": "sequencing_run",
            "exp_start_time": "2017-07-02T17:19:26Z",
            "flow_cell_id": "FAF03884",
            "heatsink_temp": "34.4726563",
            "hostname": "PC-KL-17324",
            "installation_type": "map",
            "local_firmware_file": "0",
            "msg_id": "96f773be-eb1c-473d-aa3c-b40fdad9a228",
            "operating_system": "Windows 6.2",
            "protocol_run_id": "a413384e-41eb-4264-bacc-77d95cd96685",
            "protocols_version": "1.5.12",
            "run_id": "c2d19c211888bc09d8e077df271f325c911c1010",
            "sample_id": "mt-2",
            "time_stamp": "2020-10-28T09:27:01Z",
            "usb_config": "1.1.1_ONT#MinION_fpga_1.1.0#ctrl#Auto",
            "version": "1.5.12"
        }
    }
]