
#if !defined(BASEPROCAFIELDLEVEL_H_INCLUDED)
#error "This file should only be included through BaseProcaFieldLevel.hpp"
#endif

#ifndef BASEPROCAFIELDLEVEL_IMPL_H_INCLUDED
#define BASEPROCAFIELDLEVEL_IMPL_H_INCLUDED



//do things at end of advance step, after RK4 calculation
template <class background_t, class proca_t>
void BaseProcaFieldLevel<background_t, proca_t>::specificAdvance()
{
    //check for nans
    if (m_p.nan_check){
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new,
                        EXCLUDE_GHOST_CELLS, disable_simd());
    }

};


#ifdef CH_USE_HDF5
//things to do before outputting a checkpoint file
template <class background_t, class proca_t>
void BaseProcaFieldLevel<background_t, proca_t>::prePlotLevel()
{
   
    fillAllGhosts(); 
    background_t background_init { m_p.background_params, m_dx };
    proca_t proca_field(background_init, m_p.matter_params);
    ProcaSquared<background_t> Asquared(m_dx, m_p.center, background_init);
    ChargesFluxes<proca_t, background_t> EM(background_init,m_dx, proca_field, m_p.center);

    //compute diagnostics on each cell of current level
    BoxLoops::loop(
        make_compute_pack(
            Asquared,
            EM
            ),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS
        );

};
#endif //CH_USE_HDF5

template <class background_t, class proca_t>
void BaseProcaFieldLevel<background_t, proca_t>::preTagCells()
{

}


//RHS routines used at each RK4 step
template <class background_t, class proca_t>
void BaseProcaFieldLevel<background_t, proca_t>::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                const double a_time)
{

    //Calculate right hand side with matter_t = ProcaField
    background_t background_init { m_p.background_params, m_dx };
    proca_t proca_field(background_init, m_p.matter_params);
    
    MatterEvolution<proca_t, background_t> matter_class(proca_field, background_init, m_p.sigma, m_dx, m_p.center);
    ExcisionEvolution<proca_t, background_t> excisor(m_dx, m_p.center, background_init);

    BoxLoops::loop(matter_class, a_soln, a_rhs, SKIP_GHOST_CELLS);
    BoxLoops::loop(excisor,a_soln, a_rhs, SKIP_GHOST_CELLS, disable_simd());

};




//compute tagging criteria for grid
template <class background_t, class proca_t>
void BaseProcaFieldLevel<background_t, proca_t>::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                             const FArrayBox &current_state,
                                             const FArrayBox &current_state_diagnostics)
{
    CH_TIME("BaseProcaFieldLevel::computeTaggingCriterion");

    CustomTaggingCriterion tagger(m_dx, m_level, m_p.grid_scaling*m_p.L, m_p.initial_ratio, m_p.center, m_p.extraction_params, m_p.activate_extraction, m_p.activate_ham_tagging);

    BoxLoops::loop(tagger, current_state_diagnostics, tagging_criterion, disable_simd());
}


template <class background_t, class proca_t>
void BaseProcaFieldLevel<background_t, proca_t>::specificPostTimeStep()
{
    CH_TIME("BaseProcaFieldLevel::specificPostTimeStep");

    bool first_step = (m_time == 0.); //is this the first call of posttimestep? Recall, we're calling PostTimeStep in the main function, so m_time==0 is first step

    pout() << "BaseProcaFieldLevel::specificPostTimeStep" << endl;
    int min_level = 0;
    bool at_course_timestep_on_any_level = at_level_timestep_multiple(min_level);

    //extract fluxes at specified radii
    if (m_p.activate_extraction)
    {
        
        pout() << "Extracting fluxes" << endl;
        int min_level = m_p.extraction_params.min_extraction_level();
        pout() << "min_level = " << min_level << endl;
        if (at_course_timestep_on_any_level)
        {
            pout() << "Filling ghosts" << endl;
            fillAllGhosts();
            
            pout() << "Class instantiations" << endl;
            background_t background_init { m_p.background_params, m_dx };
            proca_t proca_field(background_init, m_p.matter_params);
            ChargesFluxes<proca_t, background_t> EM(background_init,m_dx, proca_field, m_p.center);
            ExcisionDiagnostics<proca_t,background_t> excisor(background_init, m_dx, m_p.center);

            pout() << "Looping EM" << endl;
            BoxLoops::loop(
                EM,
                m_state_new, m_state_diagnostics, SKIP_GHOST_CELLS);

            pout() << "Looping excisor" << endl;
            BoxLoops::loop(
                excisor,
                m_state_diagnostics, m_state_diagnostics, SKIP_GHOST_CELLS, disable_simd()
            );

            //do extraction
            if (m_level == min_level)
            {
                CH_TIME("WeylExtraction");
                //refresh interpolator
                //fill ghosts manually
                pout() << "Refreshing interpolator" << endl;
                bool fill_ghosts = false;
                m_gr_amr.m_interpolator->refresh(fill_ghosts);
                m_gr_amr.fill_multilevel_ghosts(
                    VariableType::diagnostic, Interval(c_Edot, c_Jdot),
                    min_level);
                FluxExtraction my_extraction(m_p.extraction_params, m_dt, m_time, first_step, m_restart_time);
                pout() << "Executing query" << endl;
                my_extraction.execute_query(m_gr_amr.m_interpolator);
            }
        }
    } //end of flux block
    pout() << "Done with extraction" << endl;

    //integrate charges
    if (m_p.activate_integration && at_course_timestep_on_any_level)
    {
        pout() << "Integrating charges" << endl;
        
        //calculate densities on grid
        if ( !m_p.activate_extraction ) //did we already calculate diagnostics during extraction?
        {
            fillAllGhosts();
            pout() << "extraction not activated, so calculating charges" << endl;
            background_t background_init { m_p.background_params, m_dx };
            proca_t proca_field(background_init, m_p.matter_params);

            ChargesFluxes<proca_t, background_t> EM(background_init,m_dx, proca_field, m_p.center);
            ExcisionDiagnostics<proca_t,background_t> excisor(background_init, m_dx, m_p.center);

            BoxLoops::loop(EM,m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
            BoxLoops::loop(excisor, m_state_diagnostics, m_state_diagnostics, EXCLUDE_GHOST_CELLS, disable_simd());
        }


        if (m_level == min_level)
        {
            pout() << "We're in min_level, so perform integration" << endl;
            //setup integrator
            AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);

            //integrate densities
            double SUM_rho = amr_reductions.sum(c_rho);
            double SUM_rhoE = amr_reductions.sum(c_rhoE);
            double SUM_rhoJ = amr_reductions.sum(c_rhoJ);

            //setup output file
            SmallDataIO constraint_file(m_p.data_path + m_p.integrals_filename,
                                        m_dt, m_time, m_restart_time, SmallDataIO::APPEND, first_step
            );

            //remove duplicates
            constraint_file.remove_duplicate_time_data();

            //save to disk
            if (first_step)
            {
                constraint_file.write_header_line({"SUM_rho", "SUM_rhoE", "SUM_rhoJ"});
            }
            constraint_file.write_time_data_line({ SUM_rho, SUM_rhoE, SUM_rhoJ});
        }
         
    } //end of integration block
    pout() << "postTimeStep done" << endl;

}

#endif //BASEPROCAFIELDLEVEL_IMPL_H_INCLUDED