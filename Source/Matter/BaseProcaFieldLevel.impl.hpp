/* GRBoondi 2024
 * Please refer to LICENSE in GRBoondi's root directory.
 */

/*
    Implementation file for BaseProcaFieldLevel.hpp
*/
#if !defined(BASEPROCAFIELDLEVEL_H_INCLUDED)
#error "This file should only be included through BaseProcaFieldLevel.hpp"
#endif

#ifndef BASEPROCAFIELDLEVEL_IMPL_H_INCLUDED
#define BASEPROCAFIELDLEVEL_IMPL_H_INCLUDED

// do things at end of advance step, after RK4 calculation
template <class background_t, class proca_t>
void BaseProcaFieldLevel<background_t, proca_t>::specificAdvance()
{
    // check for nans
    if (m_p.nan_check)
    {
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new,
                       EXCLUDE_GHOST_CELLS, disable_simd());
    }
};

#ifdef CH_USE_HDF5
// things to do before outputting a checkpoint file
template <class background_t, class proca_t>
void BaseProcaFieldLevel<background_t, proca_t>::prePlotLevel()
{
    // Fill a ghosts cells. Important if we want to calculate derivatives
    fillAllGhosts();

    // initialize background and Proca classes
    background_t background_init{m_p.background_params, m_dx};
    proca_t proca_field(background_init, m_p.matter_params);

    // diagnostic class to calculated various quantities
    ProcaSquared<background_t> Asquared(m_dx, m_p.center, background_init);
    ChargesFluxes<proca_t, background_t> EM(background_init, m_dx, proca_field,
                                            m_p.center);
    DampingFieldDiagnostic z_field_diagnostic{};

    // compute diagnostics on each cell of current level
    BoxLoops::loop(make_compute_pack(Asquared, EM, z_field_diagnostic),
                   m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

    // Excise diagnostics within the excision zone.
    // Variables to excise defined by user
    const std::vector<int> vars_to_excise =
        DiagnosticVariables::convert_pairs_to_enum(
            m_p.diagnostic_excision_vars);
    ExcisionDiagnostics<proca_t, background_t> diag_excisor(
        background_init, m_dx, m_p.center, m_p.diagnostic_inner_boundary,
        m_p.diagnostic_outer_boundary, vars_to_excise);

    // excise within the excision zone
    BoxLoops::loop(diag_excisor, m_state_diagnostics, m_state_diagnostics,
                   SKIP_GHOST_CELLS, disable_simd());

    // add any other computations from the user here
    additionalPrePlotLevel();
};
#endif // CH_USE_HDF5

// RHS routines used at each RK4 step
template <class background_t, class proca_t>
void BaseProcaFieldLevel<background_t, proca_t>::specificEvalRHS(
    GRLevelData &a_soln, GRLevelData &a_rhs, const double a_time)
{

    // Calculate right hand side
    background_t background_init{m_p.background_params, m_dx};
    proca_t proca_field(background_init, m_p.matter_params);

    MatterEvolution<proca_t, background_t> matter_class(
        proca_field, background_init, m_p.sigma, m_dx, m_p.center);
    ExcisionEvolution<proca_t, background_t> excisor(
        m_dx, m_p.center, m_p.evolution_excision_width, background_init);

    BoxLoops::loop(matter_class, a_soln, a_rhs, SKIP_GHOST_CELLS);
    BoxLoops::loop(excisor, a_soln, a_rhs, SKIP_GHOST_CELLS, disable_simd());
};

// compute tagging criteria for grid
template <class background_t, class proca_t>
void BaseProcaFieldLevel<background_t, proca_t>::computeTaggingCriterion(
    FArrayBox &tagging_criterion, const FArrayBox &current_state,
    const FArrayBox &current_state_diagnostics)
{
    CH_TIME("BaseProcaFieldLevel::computeTaggingCriterion");

    // Basic tagging class
    TaggingCriterion tagger(
        m_dx, m_level, m_p.grid_scaling * m_p.L, m_p.initial_ratio, m_p.center,
        m_p.extraction_params, m_p.activate_extraction,
        m_p.activate_ham_tagging, m_p.activate_extraction_tagging);

    BoxLoops::loop(tagger, current_state_diagnostics, tagging_criterion);
}

template <class background_t, class proca_t>
void BaseProcaFieldLevel<background_t, proca_t>::specificPostTimeStep()
{
    CH_TIME("BaseProcaFieldLevel::specificPostTimeStep");

    bool first_step =
        (m_time ==
         0.); // is this the first call of posttimestep? Recall, we're calling
              // PostTimeStep in the main function, so m_time==0 is first step

    int min_level =
        m_p.extraction_params
            .min_extraction_level(); // get the minimum level for extraction, as
                                     // specified in parameter file
    bool at_course_timestep_on_any_level =
        at_level_timestep_multiple(min_level);

    // extract fluxes at specified radii
    if (m_p.activate_extraction && at_course_timestep_on_any_level &&
        m_p.num_extraction_vars > 0)
    {

        // fill the ghost cells, so we can calculate derivatives
        fillAllGhosts();

        // get enums of variables we want to extract
        const std::vector<int> vars_to_extract =
            DiagnosticVariables::convert_pairs_to_enum(m_p.extraction_vars);

        // Instantiate classes to compute the charges and fluxes
        background_t background_init{m_p.background_params, m_dx};
        proca_t proca_field(background_init, m_p.matter_params);

        // If we want to extract Edot or Jdot, we need to compute the charges
        // and fluxes on the grid
        if (DiagnosticVariables::is_variable_to_extract(c_Edot,
                                                        vars_to_extract) ||
            DiagnosticVariables::is_variable_to_extract(c_Jdot,
                                                        vars_to_extract))
        {
            ChargesFluxes<proca_t, background_t> EM(background_init, m_dx,
                                                    proca_field, m_p.center);
            // Loop over box cells and compute the charges and fluxes
            BoxLoops::loop(EM, m_state_new, m_state_diagnostics,
                           SKIP_GHOST_CELLS);
        }

        if (DiagnosticVariables::is_variable_to_extract(c_fluxLinMom,
                                                        vars_to_extract))
        {
            LinearMomConservation<proca_t, background_t> LinearMomentum(
                proca_field, background_init, m_p.linear_momentum_dir, m_dx,
                m_p.center);
            // Loop over box cells and compute the charges and fluxes
            BoxLoops::loop(LinearMomentum, m_state_new, m_state_diagnostics,
                           SKIP_GHOST_CELLS);
        }

        const std::vector<int> vars_to_excise =
            DiagnosticVariables::convert_pairs_to_enum(
                m_p.diagnostic_excision_vars);
        ExcisionDiagnostics<proca_t, background_t> excisor(
            background_init, m_dx, m_p.center, m_p.diagnostic_inner_boundary,
            m_p.diagnostic_outer_boundary, vars_to_excise);

        // excise within the excision zone
        BoxLoops::loop(excisor, m_state_diagnostics, m_state_diagnostics,
                       SKIP_GHOST_CELLS, disable_simd());

        // do extraction
        if (m_level == min_level)
        {
            CH_TIME("WeylExtraction");
            // refresh interpolator
            // fill ghosts manually
            bool fill_ghosts = false;
            m_gr_amr.m_interpolator->refresh(fill_ghosts);

            // fill multilevel ghost cells, based on the quantities we want to
            // extract, as specified in the param file
            for (auto var_enum : vars_to_extract)
            {
                m_gr_amr.fill_multilevel_ghosts(VariableType::diagnostic,
                                                Interval(var_enum, var_enum),
                                                min_level);
            };
            FluxExtraction my_extraction(
                m_flux_container, m_p.extraction_params, vars_to_extract, m_dt,
                m_time, first_step, m_restart_time, m_p.SymmetryFactor);
            my_extraction.execute_query(m_gr_amr.m_interpolator);
        }

    } // end of flux block

    // integrate charges
    if (m_p.activate_integration && at_course_timestep_on_any_level)
    {
        // get enums of variables we want to extract
        const std::vector<int> vars_to_integrate =
            DiagnosticVariables::convert_pairs_to_enum(m_p.integration_vars);

        // initialize background and matter class
        background_t background_init{m_p.background_params, m_dx};
        proca_t proca_field(background_init, m_p.matter_params);

        // calculate densities on grid
        if (!m_p.activate_extraction) // did we already fill ghosts in the
                                      // extraction block?
        {
            fillAllGhosts();
        }

        // If we want to extract any of the densities, we need to compute the
        // charges and fluxes on the grid determine if any of the variables we
        // want to integrate are calculated in ChargesFluxes class
        bool calculate_chargesfluxes = false;
        for (auto var_enum : {c_rho, c_rhoJ, c_rhoE, c_EM_squared, c_EM_trace})
        {
            calculate_chargesfluxes +=
                DiagnosticVariables::is_variable_to_extract(var_enum,
                                                            vars_to_integrate);
        }
        // if any of them are, we should compute the class over the grid
        if (calculate_chargesfluxes)
        {
            ChargesFluxes<proca_t, background_t> EM(background_init, m_dx,
                                                    proca_field, m_p.center);
            // Loop over box cells and compute the charges and fluxes
            BoxLoops::loop(EM, m_state_new, m_state_diagnostics,
                           SKIP_GHOST_CELLS);
        }

        // If we want to extract these, we need to compute the
        // LinearMomConservation class over the grid
        bool calculate_linmom = false;
        for (auto var_enum : {c_rhoLinMom, c_sourceLinMom})
        {
            calculate_linmom += DiagnosticVariables::is_variable_to_extract(
                var_enum, vars_to_integrate);
        }
        if (calculate_linmom)
        {
            LinearMomConservation<proca_t, background_t> LinearMomentum(
                proca_field, background_init, m_p.linear_momentum_dir, m_dx,
                m_p.center);
            // Loop over box cells and compute the charges and fluxes
            BoxLoops::loop(LinearMomentum, m_state_new, m_state_diagnostics,
                           SKIP_GHOST_CELLS);
        }

        const std::vector<int> vars_to_excise =
            DiagnosticVariables::convert_pairs_to_enum(
                m_p.diagnostic_excision_vars);
        ExcisionDiagnostics<proca_t, background_t> excisor(
            background_init, m_dx, m_p.center, m_p.diagnostic_inner_boundary,
            m_p.diagnostic_outer_boundary, vars_to_excise);

        // excise within the excision zone
        BoxLoops::loop(excisor, m_state_diagnostics, m_state_diagnostics,
                       SKIP_GHOST_CELLS, disable_simd());

        if (m_level == min_level)
        {
            // setup integrator
            AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);

            // setup container for solution
            std::vector<double> integrals;

            // now perform the integrals and push solution onto container
            for (auto var_enum : vars_to_integrate)
            {
                integrals.push_back(m_p.SymmetryFactor *
                                    amr_reductions.sum(var_enum));
            }

            // setup output file
            SmallDataIO constraint_file(m_p.data_path + m_p.integrals_filename,
                                        m_dt, m_time, m_restart_time,
                                        SmallDataIO::APPEND, first_step);

            // remove duplicates
            constraint_file.remove_duplicate_time_data();

            // if this is the first step, we should write a header line
            // the header name comes from DiagnosticVariables::variable_names
            if (first_step)
            {
                // first create a container for the variable names and iterate
                // over the vars_to_integrate
                std::vector<std::string> header_names;
                for (auto var_enum : vars_to_integrate)
                {
                    header_names.push_back(
                        DiagnosticVariables::variable_names[var_enum]);
                }

                // now write the header line to file
                constraint_file.write_header_line(header_names);
            }

            // now add the integrals to the file
            constraint_file.write_time_data_line(integrals);

            // update the corresponding container boject
            m_integral_container.update(m_time, integrals);
        }

    } // end of integration block

    // add any other computations from the user, via virtual function
    additionalPostTimeStep();
}

#endif // BASEPROCAFIELDLEVEL_IMPL_H_INCLUDED