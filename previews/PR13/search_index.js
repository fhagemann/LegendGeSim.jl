var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/#Modules","page":"API","title":"Modules","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Order = [:module]","category":"page"},{"location":"api/#Types-and-constants","page":"API","title":"Types and constants","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Order = [:type, :constant]","category":"page"},{"location":"api/#Functions-and-macros","page":"API","title":"Functions and macros","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Order = [:macro, :function]","category":"page"},{"location":"api/#Documentation","page":"API","title":"Documentation","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [LegendGeSim]\nOrder = [:module, :type, :constant, :macro, :function]","category":"page"},{"location":"api/#LegendGeSim.LegendGeSim","page":"API","title":"LegendGeSim.LegendGeSim","text":"LegendGeSim\n\nTemplate for Julia packages.\n\n\n\n\n\n","category":"module"},{"location":"api/#LegendGeSim.CC2","page":"API","title":"LegendGeSim.CC2","text":"CC2 circuit model (ToDo)\n\n\n\n\n\n","category":"type"},{"location":"api/#LegendGeSim.DAQ","page":"API","title":"LegendGeSim.DAQ","text":"The DAQ supertype corresponds to the component of the      DAQ and electronics setup that saves the waveform after      it triggered\n\n\n\n\n\n","category":"type"},{"location":"api/#LegendGeSim.DAQ-Tuple{PropDicts.PropDict}","page":"API","title":"LegendGeSim.DAQ","text":"DAQ(sim_conf)\n\nPropDict -> <DAQ>\n\nConstruct a DAQ supertype struct based on given simulation configuration. Type of returned instance depends on settings in <sim_conf> Currently only one type of DAQ available (GenericDAQ),     rendering this function temporarily redundant.\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.ElecChain","page":"API","title":"LegendGeSim.ElecChain","text":"The ElecChain supertype corresponds to the chain of      electronic components involved in the DAQ setup     that appear before the trigger.\n\n\n\n\n\n","category":"type"},{"location":"api/#LegendGeSim.ElecChain-Tuple{PropDicts.PropDict}","page":"API","title":"LegendGeSim.ElecChain","text":"ElecChain(sim_conf)\n\nPropDict -> <ElecChain>\n\nConstruct an ElecChain supertype struct based on given simulation configuration. Type of returned instance depends on settings in <sim_conf> Currently only one type of ElecChain available (GenericElecChain),     rendering this function temporarily redundant.\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.Environment-Tuple{PropDicts.PropDict}","page":"API","title":"LegendGeSim.Environment","text":"Simulation parameters related to the environment\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.FADC","page":"API","title":"LegendGeSim.FADC","text":"The FADC supertype corresponds to the FADC component     of the electronics chain in the DAQ setup.\n\n\n\n\n\n","category":"type"},{"location":"api/#LegendGeSim.FADC-Tuple{PropDicts.PropDict}","page":"API","title":"LegendGeSim.FADC","text":"FADC(sim_conf)\n\nPropDict -> <FADC>\n\nConstruct an FADC supertype instance based on simulation     configuration <simconf>. Type of <FADC> depends on the type specified in <simconf>     (e.g. generic, Flashcam, Struck)\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.Flashcam","page":"API","title":"LegendGeSim.Flashcam","text":"FADC with Flashcam algorithm (ToDo)\n\n\n\n\n\n","category":"type"},{"location":"api/#LegendGeSim.GenericDAQ","page":"API","title":"LegendGeSim.GenericDAQ","text":"GenericDAQ is a dummy DAQ model that stores the waveform     after it triggered, saving <baseline_length> samples before the     trigger index, and in total <nsamples> samples from start to end\n\n\n\n\n\n","category":"type"},{"location":"api/#LegendGeSim.GenericDAQ-Tuple{PropDicts.PropDict}","page":"API","title":"LegendGeSim.GenericDAQ","text":"GenericDAQ(sim_conf)\n\nPropDict -> GenericDAQ \n\nConstruct a GenericDAQ struct based on given simulation configuration \n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.GenericElecChain","page":"API","title":"LegendGeSim.GenericElecChain","text":"Generic electronics chain consisting of a preamplifier     and an FADC module\n\n\n\n\n\n","category":"type"},{"location":"api/#LegendGeSim.GenericElecChain-Tuple{PropDicts.PropDict}","page":"API","title":"LegendGeSim.GenericElecChain","text":"GenericElecChain(sim_conf)\n\nPropDict -> GenericElecChain\n\nConstruct electronics components based on simulation      configuration <sim_conf> and create a GenericChain instance     based on these components.\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.GenericFADC","page":"API","title":"LegendGeSim.GenericFADC","text":"GenericFADC is a dummy FADC model that performs sampling     of the input waveform based on the sampling interval <Δt>\n\n\n\n\n\n","category":"type"},{"location":"api/#LegendGeSim.GenericFADC-Tuple{PropDicts.PropDict}","page":"API","title":"LegendGeSim.GenericFADC","text":"GenericFADC(sim_conf)\n\nPropDict -> GenericFADC\n\nConstruct a GenericFADC instance based on simulation     configuration <sim_conf>\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.GenericPreAmp","page":"API","title":"LegendGeSim.GenericPreAmp","text":"The GenericPreAmp is a dummy PreAmp model that accounts for     effects such as decay and rise time, offset, noise and gain.\n\nThis dummy model involves no real electronics response, only gain.\n\nThis struct is currently mutable because in the case of noise modelling based on data,     the gain has to match the one in the baselines extracted from data,     so in that case the PreAmp is initialized with gain = 0, and the gain     is calculated later. I assume this is temporary, as the user should know which     gain / max_e the data was produced with and give it as simulation input.\n\n\n\n\n\n","category":"type"},{"location":"api/#LegendGeSim.GenericPreAmp-Tuple{PropDicts.PropDict}","page":"API","title":"LegendGeSim.GenericPreAmp","text":"GenericPreAmp(sim_conf)\n\nPropDict -> GenericPreAmp    \n\nConstruct a GenericPreAmp instance based on given simulation configuration <sim_conf>\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.NoiseFromData","page":"API","title":"LegendGeSim.NoiseFromData","text":"The NoiseFromData model means instead of simulating given      preamp noise, and also gain and offset, these parameters will be     accounted for or inferred from real data baselines     contained in the <baseline_catalog>\n\n\n\n\n\n","category":"type"},{"location":"api/#LegendGeSim.NoiseFromData-Tuple{PropDicts.PropDict}","page":"API","title":"LegendGeSim.NoiseFromData","text":"NoiseFromData(sim_conf)\n\nPropDict -> NoiseFromData \n\nConstruct a NoiseFromData struct based on given simulation configuration <sim_conf>\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.NoiseFromSim","page":"API","title":"LegendGeSim.NoiseFromSim","text":"The NoiseFromSim model means simulating noise from scratch     starting from fano noise of the germanium crystal     and ending with noise coming from electronics components.\n\n\n\n\n\n","category":"type"},{"location":"api/#LegendGeSim.NoiseFromSim-Tuple{PropDicts.PropDict}","page":"API","title":"LegendGeSim.NoiseFromSim","text":"NoiseFromSim(sim_conf)\n\nPropDict -> NoiseFromSim \n\nConstruct a NoiseFromSim struct based on given simulation configuration <sim_conf>\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.NoiseModel","page":"API","title":"LegendGeSim.NoiseModel","text":"The NoiseModel supertype corresponds to different methods     of simulating noise coming from the electronics chain     components.\n\nIt is not simply a string in the ElecChain struct because     it must exist as a separate instance for the stp->pss step      that is independent from the elec simulation, but has to know     what type of noise model is used to deal with fano noise.\n\n\n\n\n\n","category":"type"},{"location":"api/#LegendGeSim.NoiseModel-Tuple{PropDicts.PropDict}","page":"API","title":"LegendGeSim.NoiseModel","text":"NoiseModel(sim_config)\n\nPropDict -> <NoiseModel>\n\nConstuct a NoiseModel supertype instance based on simulation settings      given in <simconfig> Type of <NoiseModel> depends on <simconfig> settings.\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.PSSimulator","page":"API","title":"LegendGeSim.PSSimulator","text":"The PSSimulator supertype defines what type of method      is used for the simulation of electric field and weighting potential,     calculation of detector capacitance etc., as well as simulation of current pulses\n\n\n\n\n\n","category":"type"},{"location":"api/#LegendGeSim.PSSimulator-Tuple{PropDicts.PropDict}","page":"API","title":"LegendGeSim.PSSimulator","text":"PSSimulator(sim_config)\n\nPropDict -> <PSSimulator>\n\nConstruct a PSSSimulator supertype instance based on given simulation     configuration <sim_config>. Returned type depends on the simulation     method given in the config.\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.PreAmp","page":"API","title":"LegendGeSim.PreAmp","text":"The PreAmp supertype corresponds to the charge sensitive preamplifier component\n\n\n\n\n\n","category":"type"},{"location":"api/#LegendGeSim.PreAmp-Tuple{PropDicts.PropDict}","page":"API","title":"LegendGeSim.PreAmp","text":"PreAmp(sim_conf)\n\nPropDict -> <PreAmp>    \n\nConstruct a PreAmp supertype instance based on given simulation configuration <simconf>. Returned type depends on the settings in <simconf>. Currently only GenericPreAmp is available.\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.SSDSimulator","page":"API","title":"LegendGeSim.SSDSimulator","text":"Simulation method: SolidStateDetectors\n\n\n\n\n\n","category":"type"},{"location":"api/#LegendGeSim.SSDSimulator-Tuple{PropDicts.PropDict}","page":"API","title":"LegendGeSim.SSDSimulator","text":"SSDSimulator(::PropDict)\n\n-> SSDSimulator\n\nConstruct SSDSimulator instance based on simulation     configuration given in <sim_conf>.\n\nCurrently SSDSimulator does not have any parameters\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.SiggenSimulator","page":"API","title":"LegendGeSim.SiggenSimulator","text":"Simulation method: siggen\n\n\n\n\n\n","category":"type"},{"location":"api/#LegendGeSim.SiggenSimulator-Tuple{PropDicts.PropDict}","page":"API","title":"LegendGeSim.SiggenSimulator","text":"SiggeSimulator(sim_conf)\n\nPropDict -> SiggenSimulator\n\nConstruct SiggeSimulator instance based on simulation     configuration given in <sim_conf>.\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.Struck","page":"API","title":"LegendGeSim.Struck","text":"FADC with Struck algorithm (ToDo)\n\n\n\n\n\n","category":"type"},{"location":"api/#LegendGeSim.TrapFilter","page":"API","title":"LegendGeSim.TrapFilter","text":"The TrapFilter trigger is a dummy simulation of the trapezoidal filter     algorithm to produce a trigger if the waveform amplitude passes     a certain threshold.\n\nThis struct is currently mutable because in case of NoiseFromData modelling,     we need to calculate the threshold post-factum after analysing the noise     levels in data.    \n\n\n\n\n\n","category":"type"},{"location":"api/#LegendGeSim.TrapFilter-Tuple{PropDicts.PropDict}","page":"API","title":"LegendGeSim.TrapFilter","text":"TrapFilter(sim_conf)\n\nPropDict -> TrapFilter \n\nConstruct a TrapFilter instance based on simulation configuration given in <sim_conf>.\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.Trigger","page":"API","title":"LegendGeSim.Trigger","text":"The Trigger supertype corresponds to the trigger component in DAQ.\n\n\n\n\n\n","category":"type"},{"location":"api/#LegendGeSim.Trigger-Tuple{PropDicts.PropDict}","page":"API","title":"LegendGeSim.Trigger","text":"Trigger(sim_conf)\n\nPropDict -> <Trigger>    \n\nConstruct a Trigger supertype instance based on settings given in <sim_conf>. The returned type depends on the given settings.\n\nCurrently only TrapFilter type is implemented.\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.add_tail_and_baseline-Tuple{RadiationDetectorSignals.RDWaveform, LegendGeSim.GenericFADC, LegendGeSim.GenericDAQ}","page":"API","title":"LegendGeSim.add_tail_and_baseline","text":"add_tail_and_baseline(wf, fadc, daq)\n\nRDWaveform, GenericFADC, GenericDAQ -> RDWaveform  \n\nExtend tail and add baseline of <wf> based on <fadc> and <daq> parameters.\n\nExtended conservatively to ensure enough samples for the future     sliding trap filter window and the resulting DAQ baseline.\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.baseline_catalog-Tuple{AbstractString}","page":"API","title":"LegendGeSim.baseline_catalog","text":"baseline_catalog(raw_filename)\n\nAbstractString -> Table \n\nLook up stored table of baselines corresponding to given raw data <raw_filename>. If does not exist, construct such table.\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.baseline_catalog-Tuple{TypedTables.Table}","page":"API","title":"LegendGeSim.baseline_catalog","text":"baseline_catalog(raw_table)\n\nTable -> Table \n\nConstruct table of baselines extracted from the waveforms     contained in the given raw data table <raw_table>\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.construct_ssd_simulation-Tuple{PropDicts.PropDict, LegendGeSim.Environment, LegendGeSim.SSDSimulator}","page":"API","title":"LegendGeSim.construct_ssd_simulation","text":"construct_ssd_simulation(det_meta::AbstractString, env::Environment, sim_settings::SSDSimulator)\n\nConstruct a SolidStateDetectors.Simulation based on geometry as given in LEGEND metadata det_meta and on envorinmental settings specified in env and on simulational settings specified in sim_settings.\n\nToDo: Actually use env. \n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.differentiate-Tuple{RadiationDetectorSignals.RDWaveform}","page":"API","title":"LegendGeSim.differentiate","text":"differentiate(wf)\n\nRDWaveform -> RDWaveform    \n\nDifferentiate a waveform using a Biquad filter\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.extend_baseline-Tuple{RadiationDetectorSignals.RDWaveform, RadiationDetectorSignals.RDWaveform}","page":"API","title":"LegendGeSim.extend_baseline","text":"extend_baseline(baseline, wf)\n\nRDWaveform, RDWaveform -> RDWaveform\n\nTake a given <baseline> and extend it to match the length (and sampling)             of the given waveform <wf>\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.fano_noise-Tuple{TypedTables.Table, PropDicts.PropDict, LegendGeSim.Environment, LegendGeSim.NoiseFromData}","page":"API","title":"LegendGeSim.fano_noise","text":"fano_noise(events, ::PropDict, ::Environment, ::NoiseFromData)\n\nTable -> Table    \n\nDo nothing since we do not need to simulate fano noise separately when      using data baselines to account for noise levels.\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.fano_noise-Tuple{TypedTables.Table, PropDicts.PropDict, LegendGeSim.Environment, LegendGeSim.NoiseFromSim}","page":"API","title":"LegendGeSim.fano_noise","text":"fano_noise(events, det_meta, env, ::NoiseFromSim)\n\nTable, PropDict, Environment -> Table \n\nCalculate fano noise level based on the detector specification provided in     LEGEND metadata <det_meta> and environment settings provided in <env>,     and add it to given <events>\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.filename-Tuple{Any}","page":"API","title":"LegendGeSim.filename","text":"filename(path)\n\nSring -> String\n\nExtract core name of the file from path \n\nE.g. filename(\"/some/path/to/filename.ext\") -> \"filename\"\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.group_by_column-Tuple{TypedTables.Table, Symbol}","page":"API","title":"LegendGeSim.group_by_column","text":"group_by_column(table, colname)\n\nTypedTables.Table, Symbol -> TypedTables.Table \n\nGroup table <table> by given column <colname>.\n\nI think this function already exists somewhere. Gotta ask Oliver and maybe use that package instead.\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.load_config-Tuple{AbstractString, AbstractString}","page":"API","title":"LegendGeSim.load_config","text":"load_config(input_file, det_metadata, sim_config_file)\n\nAbstractString or Table, AbstractString, AbstractString -> PropDict\n\nMerge simulation inputs and settings into one simulation config to pass around\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.meta2siggen-Tuple{PropDicts.PropDict, LegendGeSim.Environment}","page":"API","title":"LegendGeSim.meta2siggen","text":"meta2siggen(meta, env)\n\nPropDict, Environment -> Dict    \n\nConstruct siggen type config in a form of Dict based on     geometry as given in LEGEND metadata <meta>     and environment settings given in <env>.\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.pet_to_raw-Tuple{PropDicts.PropDict}","page":"API","title":"LegendGeSim.pet_to_raw","text":"pet_to_raw(pet_file_path, sim_config, config_name)\n\nAbstractString, PropDict, AbstractString -> Table\n\n[WIP] Full simulation chain pet->stp->pss->raw     based on simulated energy depositions contained in the HDF5 file     found in <petfilepath> and simulation settings given in <sim_config>.\n\nThe output name <config_name> is used to construct filenames for cached     simulation files (currently the same as the simulation config basename)\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.pet_to_stp-Union{Tuple{S}, Tuple{TypedTables.Table, S}} where S","page":"API","title":"LegendGeSim.pet_to_stp","text":"pet_to_stp(pet_table, detector_SSD)\n\nTable, SolidStateDetectors.Simulation -> Table\n\nConstruct a table with \"stepping info\" based on     position-energy-time information of simulated energy depositions     given in <pettable> and detector geometry given in <detectorSSD>.\n\nCurrent steps include:     - removing events outside of detector PV     - clustering     - grouping by detector     - removing events reconstructed to be outisde of the detector\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.preamp_gain-Tuple{LegendGeSim.GenericPreAmp, LegendGeSim.NoiseFromSim}","page":"API","title":"LegendGeSim.preamp_gain","text":"preamp_gain(preamp, ::NoiseFromSim)\n\nGenericPreAmp -> Float64\n\nDo nothing and return original value of gain parameter in <preamp>,     since when NoiseFromSim model is used, gain is provided by the user     (or calculated based on user given max energy)    \n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.preamp_gain-Tuple{LegendGeSim.PreAmp, LegendGeSim.NoiseFromData}","page":"API","title":"LegendGeSim.preamp_gain","text":"preamp_gain(preamp, noise_model)\n\nGenericPreAmp, NoiseFromData -> Real\n\nCalculate gain of <preamp> based on its max energy and      the offset observed in data baselines contained in <noise_model>\n\nIn principle we should not do this in this simulation, and user has to provide      the precise parameters of the electronics chain used in producing data      the baselines from which are contained in <noise_model>.\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.process_waveforms-Tuple{TypedTables.Table, LegendGeSim.PSSimulator, LegendGeSim.GenericElecChain, LegendGeSim.Trigger, LegendGeSim.DAQ, LegendGeSim.NoiseModel}","page":"API","title":"LegendGeSim.process_waveforms","text":"process_waveforms(pss_table, elec_chain, trigger, daq, noise_model)\n\nTable, GenericElecChain, Trigger, DAQ, NoiseModel -> Table  \n\nSimulate effects of the electronics chain <elecchain>, <trigger> and <daq> settings     on the waveforms contained in <psstable>. Noise is simulated based on the given <noise_model>. The output table contains the resulting DAQ waveforms, simulated online energy after tigger,     baselines and their RMS, needed for the raw tier table.\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.propdict-Tuple{AbstractString}","page":"API","title":"LegendGeSim.propdict","text":"propdict(json_file)\n\nAbstractString -> PropDict \n\nConstruct a PropDict based on given <json_file>.\n\nWhy can't I name it function PropDict()?\nWhy doesn't this already exist in PropDicts?\n\nI find the PropDicts.read(PropDict, String) format kinda cumbersome\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.pss_to_raw-Tuple{TypedTables.Table, TypedTables.Table, LegendGeSim.PSSimulator, LegendGeSim.ElecChain, LegendGeSim.Trigger, LegendGeSim.DAQ, LegendGeSim.NoiseModel}","page":"API","title":"LegendGeSim.pss_to_raw","text":"pss_to_raw(pss_table, pss_truth, simulation_settings, elec_chain, trigger, daq, noise_model)\n\nTable, Table, SSDSimulator, ElecChain, Trigger, DAQ, NoiseModel -> Table\n\nSimulate effects of the electronics chain <elecchain>, <trigger> and <daq> settings     on the waveforms contained in <psstable>. Noise is simulated based on the given <noisemodel>. Construct a table with the format identical to data raw tier. Currently timing information in <psstruth> is used for dummy timestamps in the output     raw tier table.\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.read_pet-Tuple{AbstractString}","page":"API","title":"LegendGeSim.read_pet","text":"AbstractString -> Table\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.remove_negative-Tuple{Real}","page":"API","title":"LegendGeSim.remove_negative","text":"remove_negative(value)\n\nReal -> Real \n\nReplace negative values by zeros.\n\nUsed to remove a glitch in SSD pulses as a quick fix     while the glitch is being debugged.\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.siggen_config-Tuple{PropDicts.PropDict, LegendGeSim.Environment, LegendGeSim.SiggenSimulator, AbstractString}","page":"API","title":"LegendGeSim.siggen_config","text":"siggen_config(meta, env, siggen_sim, config_name)\n\nPropDict, Environment, SiggenSimulator, String -> String, String    \n\nConstruct Siggen/Fieldgen config based on     geometry as given in LEGEND metadata <meta>,     environment settings given in <env>,     and fieldgen settings contained in <siggen_sim>.\n\nThe resulting siggen config will be cached based on given <config_name> The function returns the path to the cached siggen config file, and     the relative path to the fieldgen WP file (to check if already exists later)\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.simulate-Tuple{RadiationDetectorSignals.RDWaveform, Int64, LegendGeSim.GenericDAQ}","page":"API","title":"LegendGeSim.simulate","text":"simulate(wf, trigger_index, daq)\n\nRDWaveform, Int, DAQ -> RDWaveform\n\nSimulate how the DAQ stores the waveform after it receives a trigger.\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.simulate-Tuple{RadiationDetectorSignals.RDWaveform, LegendGeSim.ElecChain, LegendGeSim.NoiseFromSim}","page":"API","title":"LegendGeSim.simulate","text":"simulate(wf, elec_chain, ::NoiseFromSim)\n\nRDWaveform, ElecChain -> RDWaveform\n\nSimulate effects of the electronics chain on the waveform     modelling noise based on each component.\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.simulate-Tuple{RadiationDetectorSignals.RDWaveform, LegendGeSim.ElecChain}","page":"API","title":"LegendGeSim.simulate","text":"simulate(wf, elec_chain)\n\nRDWaveform, ElecChain -> RDWaveform\n\nSimulate effects of the electronics chain components on the waveform.\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.simulate-Tuple{RadiationDetectorSignals.RDWaveform, LegendGeSim.GenericElecChain, LegendGeSim.NoiseFromData}","page":"API","title":"LegendGeSim.simulate","text":"simulate(wf, elec_chain, noise_model)\n\nRDWaveform, GenericElecChain, NoiseFromData -> RDWaveform\n\nSimulate effects of the electronics chain on the waveform     modelling noise based on baselines extracted from data     (note: also takes preamp offset into account)\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.simulate-Tuple{RadiationDetectorSignals.RDWaveform, LegendGeSim.GenericFADC}","page":"API","title":"LegendGeSim.simulate","text":"simulate(wf, fadc)\n\nRDWaveform, GenericFADC -> RDWaveform\n\nSimulate effects of FADC module <fadc> on the waveform <wf>     (sampling and digitization)\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.simulate-Tuple{RadiationDetectorSignals.RDWaveform, LegendGeSim.GenericPreAmp}","page":"API","title":"LegendGeSim.simulate","text":"simulate(wf, preamp)\n\nRDWaveform, GenericPreAmp -> RDWaveform\n\nSimulate the effecs of the preamp on the waveform.\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.simulate-Tuple{RadiationDetectorSignals.RDWaveform, LegendGeSim.TrapFilter}","page":"API","title":"LegendGeSim.simulate","text":"simulate(wf, trap_filter)\n\nRDWaveform, TrapFilter -> Int, Real \n\nSimulate the results of applying <trapfilter> to <wf> Returns the index on which <wf> triggered, and the estimated     online energy based on <trapfilter> output.\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.simulate_detector-Tuple{PropDicts.PropDict, LegendGeSim.Environment, AbstractString, LegendGeSim.SiggenSimulator}","page":"API","title":"LegendGeSim.simulate_detector","text":"simulate_detector(det_meta, env, config_name, ps_simulator)\n\nPropDict, Environment, AbstractString, Siggen -> SigGenSetup       \n\nConstruct a fieldgen/siggen configuration file     according to geometry as given in LEGEND metadata <detmeta>     and environment settings given in <env>. Look up fieldgen generated electric field and weighting potential files     based on given <configname>, and if not found, call fieldgen.\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.simulate_noise-Tuple{RadiationDetectorSignals.RDWaveform, LegendGeSim.NoiseFromData}","page":"API","title":"LegendGeSim.simulate_noise","text":"simulate_noise(wf, noise_model)\n\nRDWaveform, NoiseFromData -> RDWaveform\n\nSimulate noise and offset by picking a random baseline from the baseline catalog     contained in <noise_model> and slapping it on top of the given waveform <wf>.\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.simulate_noise-Tuple{RadiationDetectorSignals.RDWaveform, LegendGeSim.PreAmp}","page":"API","title":"LegendGeSim.simulate_noise","text":"simulate_noise(wf, preamp)\n\nRDWaveform, PreAmp -> RDWaveform\n\nSimulate effects of the preamplifier <preamp> on the given waveform <wf>.\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.simulate_signal-Tuple{AbstractVector, AbstractVector, MJDSigGen.SigGenSetup}","page":"API","title":"LegendGeSim.simulate_signal","text":"simulate_wf(pos, edep, siggen_setup)\n\nAbstractVector, AbstractVector, SigGenSetup -> Vector{Float32}\n\nSimulate a signal from energy depositions <edep> with corresponding positions <pos>     based on a given <siggen_setup>    \n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.simulate_waveforms-Tuple{TypedTables.Table, MJDSigGen.SigGenSetup}","page":"API","title":"LegendGeSim.simulate_waveforms","text":"simulate_waveforms(stp_events, detector)\n\nTable, AbstractString -> Table, Table     \n\nSimulate pulses based on events given in <stp_events>     using Siggen with settings given in <detector>\n\nConstructs and returns a table with resulting pulses and a pss truth table     (may be abolished in the future as unnecessary)    \n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.simulate_waveforms-Tuple{TypedTables.Table, SolidStateDetectors.Simulation}","page":"API","title":"LegendGeSim.simulate_waveforms","text":"simulate_waveforms(stp_events, detector)\n\nTable, SolidStateDetectors.Simulation -> Table, Table     \n\nSimulate pulses based on events given in <stp_events>     using the given SSD detector simulation instance <detector>\n\nConstructs and returns a table with resulting pulses and a pss truth table     (may be abolished in the future as unnecessary)    \n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.stp_to_pss-Tuple{TypedTables.Table, PropDicts.PropDict, LegendGeSim.Environment, LegendGeSim.SSDSimulator, LegendGeSim.NoiseModel, AbstractString}","page":"API","title":"LegendGeSim.stp_to_pss","text":"stp_to_pss(stp_table, det_meta, env, ps_simulator, noise_model, config_name)\n\nTable, PropDict, Env, PSSimulator, NoiseModel, AbstractString -> Table, Table    \n\nSimulate waveforms based on stepping info given in <stptable>,      LEGEND detector metadata <detmeta>, environment settings given in <env>,     pulse shape simulation method <pssimulator>, and <noisemodel>\n\nThe output is a table with simulated pulses, and a table with simulation truth     (may be abolished in the future since basically corresponds to stp table)\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.trap_filter-Tuple{AbstractVector, Int64, LegendGeSim.TrapFilter}","page":"API","title":"LegendGeSim.trap_filter","text":"trap_filter(wf_values, sample_idx, trap_filter)\n\nAbstractVector, Int, TrapFilter -> Real, Bool \n\nSimulate the output of the trapezoidal filter algorithm     applied to <wfvalues> when the first sliding window starts     at <sampleidx>. The trapezoidal algorithm parameters are contained in <trap_filter. The returned values are the online energy estimate corresponding to the     difference between the 1st and the 3rd sliding window, and a Bool     correspondong to whether there is a trigger issued or not.    \n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.trigger_threshold-Tuple{LegendGeSim.Trigger, LegendGeSim.GenericPreAmp, LegendGeSim.NoiseFromData}","page":"API","title":"LegendGeSim.trigger_threshold","text":"trigger_threshold(trigger, preamp, noise_model)\n\nTrigger, GenericPreAmp, NoiseFromData -> Real \n\nIn NoiseFromData setting, if trigger threshold in keV is not provided,     it is calculated based on noise levels in the baselines contained in     <noise_model>. Otherwise the final threshold is calculated based on <preamp> gain.\n\n\n\n\n\n","category":"method"},{"location":"api/#LegendGeSim.trigger_threshold-Tuple{LegendGeSim.Trigger, LegendGeSim.GenericPreAmp, LegendGeSim.NoiseFromSim}","page":"API","title":"LegendGeSim.trigger_threshold","text":"trigger_threshold(trigger, preamp, ::NoiseFromSim)\n\nTrigger, GenericPreAmp -> Real \n\nIn NoiseFromSim setting, non-zero trigger threshold in keV MUST be given by the user. Calculate final threshold based the value contained in <trigger>, and <preamp> gain.\n\n\n\n\n\n","category":"method"},{"location":"LICENSE/#LICENSE","page":"LICENSE","title":"LICENSE","text":"","category":"section"},{"location":"LICENSE/","page":"LICENSE","title":"LICENSE","text":"using Markdown\nMarkdown.parse_file(joinpath(@__DIR__, \"..\", \"..\", \"LICENSE.md\"))","category":"page"},{"location":"#LegendGeSim.jl","page":"Home","title":"LegendGeSim.jl","text":"","category":"section"}]
}
