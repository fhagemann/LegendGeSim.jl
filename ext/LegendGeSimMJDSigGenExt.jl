module LegendGeSimMJDSigGenExt

    using MJDSigGen

    using PropDicts
    using SolidStateDetectors


    include("mjdsiggen_utils.jl")


    ####################################
    ### FieldGen
    ####################################

    """
        simulate_detector(det_meta, env, config_name, ps_simulator)
        
    PropDict, Environment, AbstractString, Siggen -> SigGenSetup       

    Construct a fieldgen/siggen configuration file
        according to geometry as given in LEGEND metadata <det_meta>
        and environment settings given in <env>.
    Look up fieldgen generated electric field and weighting potential files
        based on given <config_name>, and if not found, call fieldgen.
    """
    function LegendGeSim.simulate_detector(det_meta::PropDict, env::Environment, simulator::SiggenSimulator;
            overwrite::Bool = false)

        @info "Extension"

        # if no cached name given, force simulation from scratch (or else might read past "tmp" file)
        if simulator.cached_name == ""
            overwrite = true
        end

        # returns the name of the resulting siggen config file
        # and the name of (already or to be) generated weighting potential file
        # ToDo: don't create if already present already at this stage -> check in siggen_config() if name exists and just return name
        siggen_config_name, fieldgen_wp_name = siggen_config(det_meta, env, simulator; overwrite)
        fieldgen_wp_name = joinpath("cache", fieldgen_wp_name)

        # if the WP file with such a name does not exist yet or want to overwrite it...
        if !isfile(fieldgen_wp_name) || overwrite
            imp_filename, offset =
                if simulator.crystal_metadata_path != ""
                    # create impurity input on the fly
                    impurity_density_model(det_meta, simulator.crystal_metadata_path, simulator)
                else
                    # user provided .dat file and corresponding offset
                    # (or if nothing is provided this will be "" and -1)
                    simulator.impurity_profile, simulator.offset_in_mm
                end
            #...call fieldgen -> will write the wp file
            @info "_|~|_|~|_|~|_ Fieldgen simulation"
            fieldgen(siggen_config_name; impurity_profile=imp_filename, offset_mm=offset)
            @info "_|~|_|~|_|~|_ Fieldgen simulation complete"
        else
            #...do nothing, siggen will later read the files based on the generated conifg
            @info "Cached simulation found. Reading cached fields from $fieldgen_wp_name"
        end

        # a SigGenSetup object
        SigGenSetup(siggen_config_name)
    end

end

"""
    simulate_waveforms(stp_events, detector)

Table, AbstractString -> Table, Table     

Simulate pulses based on events given in <stp_events>
    using Siggen with settings given in <detector>

Constructs and returns a table with resulting pulses and a pss truth table
    (may be abolished in the future as unnecessary)    
"""
function LegendGeSim.simulate_waveforms(stp_events::Table, detector::SigGenSetup, ::SiggenSimulator)
    T = Float32 # This should be somehow defined and be passed properly
    @info("~.~.~ Siggen")

    nevents = length(stp_events)
    wf_array = Array{RDWaveform}(undef, nevents)
    for i in 1:nevents
        # println("$i/$nevents")
        if(i % 100  == 0) println("$i/$nevents") end
        signal::Vector{Float32} = simulate_signal(stp_events[i].pos, stp_events[i].edep, detector)
        time = range(T(0)u"ns", step = detector.step_time_out*u"ns", length = length(signal))

        # push!(waveforms, signal)
        wf_array[i] = RDWaveform(time, signal)
    end

    ## Oliver said this should be the more efficient way of doing it
    ## as opposed to an array of separate RDWaveform instances
    # Δt = setup.step_time_out*u"ns"
    # t_end = length(waveforms[1]) * Δt
    # times = fill(0u"ns" : Δt : t_end, length(waveforms))
    # ArrayOfRDWaveforms((times, VectorOfSimilarVectors(waveforms)))

    waveforms = ArrayOfRDWaveforms(wf_array)
    # why am I doing this?
    waveforms = ArrayOfRDWaveforms((waveforms.time, VectorOfSimilarVectors(waveforms.signal)))    

    pss_table = Table(
        channel = [1 for idx in 1:length(waveforms)], # lists of ADCs that triggered, 1 for HADES all the time
        ievt = stp_events.evtno,
        waveform = waveforms
    )

    # using SSD we get mc truth from SSD output
    # with siggen, let's just return the input
    # maybe mc truth will not be propagated to (sim)raw anyway
    pss_truth = Table(
        detno = stp_events.detno,
        edep = stp_events.edep,
        ievt = stp_events.evtno,
        pos = stp_events.pos,
        thit = stp_events.thit        
    )

    pss_table, pss_truth



    """
        simulate_wf(pos, edep, siggen_setup)

    AbstractVector, AbstractVector, SigGenSetup -> Vector{Float32}

    Simulate a signal from energy depositions <edep> with corresponding positions <pos>
        based on a given <siggen_setup>    
    """
    function LegendGeSim.simulate_signal(pos::AbstractVector, edep::AbstractVector, siggen_setup::SigGenSetup)
        # this is not so nice, think about it
        signal = zeros(Float32, siggen_setup.ntsteps_out)

        # round to siggen crystal grid precision
        # a = string(siggen_setup.xtal_grid) # crystal grid precision e.g. 0.1 (mm)
        # ndig = length(a[findfirst('.', a)+1:end]) # number of digits after .

        for i in 1:length(pos)
            # pos_rounded = round.(ustrip.(pos[i]), digits=ndig) # strip of units and round to siggen precision
            # in SSD the output is always in eV -> convert to eV
            signal = signal .+ MJDSigGen.get_signal!(siggen_setup, Tuple(ustrip.(pos[i]))) * ustrip(uconvert(u"eV", edep[i]))
            # signal = signal .+ MJDSigGen.get_signal!(siggen_setup, Tuple(pos_rounded)) * ustrip(uconvert(u"eV", edep[i]))

            # somehow this doesn't work
            # MJDSigGen.get_signal!(signal, siggen_setup, Tuple(ustrip.(pos[i]))) * ustrip(uconvert(u"eV", edep[i]))
        end

        signal
    end

    function LegendGeSim.capacitance_matrix(sim::SigGenSetup) 
        c = sim.capacitance * u"pF"
        [c -c; -c c]
    end

end



