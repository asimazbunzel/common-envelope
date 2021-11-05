! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use chem_def
      
      implicit none

      include 'test_suite_extras_def.inc'

      integer :: redo_count

      logical :: dbg = .false.

      ! lxtra(lx_found_zams) is true if star reaches ZAMS according to Aaron Dotter's definition
      integer, parameter :: lx_found_zams = 1
      ! lxtra(lx_have_convective_envelope) is true when convective envelope mass is above 20% in mass
      integer, parameter :: lx_have_convective_envelope = 2
      ! lxtra(lx_hydro_on) is true when hydrodynamics is on
      integer, parameter :: lx_hydro_on = 3
      ! lxtra(lx_ce_on) is true when have_convective_envelope and after n * tdyn of envelope
      integer, parameter :: lx_ce_on = 4

      ! xtra(x_X_central_initial) is the value of [H1] at start of runs
      integer, parameter :: x_X_central_initial = 1
      ! xtra(x_dyn_time) is the dynamical timescale of the envelope in sec
      integer, parameter :: x_dyn_time = 2
      ! xtra(x_time_start_ce) is the time at the start of CE in sec
      integer, parameter :: x_time_start_ce = 3
      ! xtra(x_r_acc) is the position of the accretor measured from the center of the donor in cm
      integer, parameter :: x_r_acc = 4
      ! xtra(x_dr) is the change in position of the accretor in a single timestep in cm
      integer, parameter :: x_dr = 5
      ! xtra(x_bondi_radius) is the bondi accretion radius in cm
      integer, parameter :: x_bondi_radius = 6
      ! xtra(x_de) is the energy deposited in a single timestep in erg
      integer, parameter :: x_de = 7
      ! xtra(x_v_circ) is the circular velocity of the accretor in cm/sec
      integer, parameter :: x_v_circ = 8
      ! xtra(x_cs) is the sound speed at the position of the accretor
      integer, parameter :: x_cs = 9
      ! xtra(x_penetration_coeff) is the coeff of penetration, i.e., accretor + bondi inside donor
      integer, parameter :: x_penetration_coeff = 10
      ! xtra(x_eps_rho) is the density scale height
      integer, parameter :: x_eps_rho = 11
      ! xtra(x_v_esc) is the escape velocity reached inside the star
      integer, parameter :: x_v_esc = 12
      ! xtra(x_r_at_v_esc) is the position where escape velocity is reached inside the star
      integer, parameter :: x_r_at_v_esc = 13
      ! xtra(x_v_at_v_esc) is the velocity where escape velocity is reached inside the star
      integer, parameter :: x_v_at_v_esc = 14
      ! xtra(x_initial_energy) is the envelope energy available at start of ce
      integer, parameter :: x_initial_energy = 30


      real(dp) :: m_acc, r_acc, fraction_of_r_donor
      real(dp) :: energy_multiplier, cs_multiplier
      logical :: use_modified_HLA, remove_surface
      integer :: num_dyn

      ! these routines are called by the standard run_star check_model
      contains
      
      include 'test_suite_extras.inc'

      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         m_acc = s% x_ctrl(1) * Msun  ! in g
         fraction_of_r_donor = s% x_ctrl(2)
         energy_multiplier = s% x_ctrl(3)
         cs_multiplier = s% x_ctrl(4)
         use_modified_HLA = s% x_logical_ctrl(1)
         remove_surface = s% x_logical_ctrl(2)
         num_dyn = s% x_integer_ctrl(1)
         if (dbg) then
            write(*,'(a)')
            write(*,'(a40, f26.16)') 'accretor mass', m_acc / Msun
            write(*,'(a40, f26.16)') 'fraction of r_donor', fraction_of_r_donor
            write(*,'(a)')
         end if

         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
         s% other_wind => brott_wind
         s% other_energy => edep
         s% other_before_struct_burn_mix => before_struct_burn_mix_for_edep

         ! the extras functions in this file will not be called
         ! unless you set their function pointers as done below.
         ! otherwise we use a null_ version which does nothing (except warn).
         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns
         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

      end subroutine extras_controls

      
      include 'winds_utils.inc'


      ! extra_heat units: erg/cm/sec
      subroutine edep(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type(star_info), pointer :: s
         integer :: nz
         integer :: k, k_acc, k_bot, k_top
         real(dp) :: v_circ, cs, bondi_radius, mach_number
         real(dp) :: dm_bondi, mean_rho
         real(dp) :: penetration_coeff
         real(dp) :: cross_section_for_drag
         real(dp) :: eps_rho, f_m
         real(dp) :: de, dr

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         dr = 0d0
         de = 0d0
         nz = s% nz
         s% extra_heat(:) = 0d0

         ! calculate sound speed at the position of the wd
         k_acc = 1
         cs = s% csound(1)
         do k = 1, nz
            if (s% r(k) <= r_acc) then
               k_acc = k
               cs = s% csound(k_acc)
               exit
            end if
         end do

         ! relative orbital velocity of WD assuming circular orbit
         v_circ = orbital_velocity(s% m(k_acc) + m_acc, r_acc)

         ! eval bondi accretion radius
         bondi_radius = bondi_accretion_radius(m_acc, cs, v_circ)

         ! find bottom & top cell numbers that matches influence area of bondi accretion radius
         call get_influence_cell_k(nz, s% r, r_acc, bondi_radius, k_bot, k_top)
         if (dbg) then
            write(*,'(a40, i14)') 'model_number', s% model_number
            write(*,'(a40, i14)') 'k_acc', k_acc
            write(*,'(a40, i14)') 'k_bot', k_bot
            write(*,'(a40, i14)') 'k_top', k_top
            write(*,'(a40, f26.16)') 'R_bot', s% r(k_bot) / Rsun
            write(*,'(a40, f26.16)') 'R_top', s% r(k_top) / Rsun
            write(*,'(a40, f26.16)') 'bondi_radius', bondi_radius / Rsun
         end if

         ! mass inside this region to get mean density
         dm_bondi = sum(s% dm(k_top:k_bot))
         mean_rho = dot_product(s% rho(k_top:k_bot), s% dm(k_top:k_bot)) / dm_bondi
         if (dbg) then
            write(*,'(a40, f26.16)') 'dm_bondi', dm_bondi / Msun
            write(*,'(a40, f26.16)') 'mean_rho', mean_rho
         end if

         ! before evaluating the cross-section, check if wd + bondi_radius are completely inside
         ! the donor star, else calculate the penetration coeff.
         penetration_coeff = eval_penetration(bondi_radius, s% r(1), r_acc)
         if (dbg) write(*,'(a40, f26.16)') 'penetration_coeff', penetration_coeff / Rsun

         ! cross-section for drag, depending on penetration_coeff
         if (penetration_coeff >= 0d0 .and. (r_acc + bondi_radius) > s% r(1)) then
            cross_section_for_drag = engulfed_area(penetration_coeff, bondi_radius)
         else
            cross_section_for_drag = pi * pow2(bondi_radius)
         end if
         if (dbg) write(*,'(a40, f26.16)') 'cross-section', cross_section_for_drag / pow2(Rsun)

         ! get energy to deposit from grav drag which depends on the mach number
         ! (i.e., on the velocity of the WD wrt the sound speed)
         call eval_drag_force(s% m(k_acc), m_acc, cross_section_for_drag, mean_rho, &
            s% dt, r_acc, de, dr)

         mach_number = v_circ / cs
         if (dbg) then
            write(*,'(a40, 1pd26.16)') 'orbital velocity', v_circ / 1d5
            write(*,'(a40, 1pd26.16)') 'sound speed', cs / 1d5
            write(*,'(a40, f26.16)') 'mach number', mach_number
         end if

         ! check if use modified HLA, or classic
         if (use_modified_HLA) then
            call macleod_fit(s, k_acc, bondi_radius, eps_rho, f_m)
         else
            eps_rho = -99d0
            f_m = 1d0
         end if

         if (dbg) write(*,'(a40, f26.16)') 'f_m', f_m

         ! update dr and dt for density inhomogeneities
         dr = dr * f_m
         de = de * f_m

         ! inject energy
         do k = 1, nz
            s% extra_heat(k) = energy_multiplier * de * exp(-pow2((s% r(k)-r_acc) / bondi_radius)) &
               / dm_bondi / s% dt
         end do

         if (dbg) then
            write(*,'(a40, 1pd26.16)') 'de', de
            write(*,'(a40, 1pd26.16)') 'edep', sum(s% extra_heat(1:nz))
            write(*,'(a40, f26.16)') 'dr', dr / Rsun
            write(*,'(a40, f26.16)') 'R_phot', s% photosphere_r
            write(*,'(a40, f26.16)') 'R_surface', s% r(1) / Rsun
            write(*,'(a40, f26.16, /)') 'r_acc', r_acc / Rsun
         end if

         ! store some info for output
         s% xtra(x_r_acc) = r_acc
         s% xtra(x_dr) = dr
         s% xtra(x_bondi_radius) = bondi_radius
         s% xtra(x_de) = de
         s% xtra(x_v_circ) = v_circ
         s% xtra(x_cs) = cs
         s% xtra(x_penetration_coeff) = penetration_coeff
         s% xtra(x_eps_rho) = eps_rho

      end subroutine edep


      include 'edep_utils.inc'


      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call test_suite_startup(s, restart, ierr)

         if (.not. restart) then
            call alloc_extra_info(s)
            s% lxtra(lx_found_zams) = .false.
            s% lxtra(lx_have_convective_envelope) = .false.
            s% lxtra(lx_hydro_on) = .false.
            s% lxtra(lx_ce_on) = .false.
            s% xtra(x_X_central_initial) = s% center_h1
            s% xtra(x_dyn_time) = 0d0
            s% xtra(x_time_start_ce) = 0d0
            s% xtra(x_r_acc) = 0d0
            s% xtra(x_dr) = 0d0
            s% xtra(x_bondi_radius) = 0d0
            s% xtra(x_de) = 0d0
            s% xtra(x_v_circ) = 0d0
            s% xtra(x_cs) = 0d0
            s% xtra(x_penetration_coeff) = 0d0
            s% xtra(x_eps_rho) = 0d0
            s% xtra(x_v_esc) = 0d0
            s% xtra(x_r_at_v_esc) = 0d0
            s% xtra(x_v_at_v_esc) = 0d0
         else
            call unpack_extra_info(s)
            ! if restart past ZAMS
            if (s% lxtra(lx_found_zams)) then
               s% job% pgstar_flag = .true.
               s% use_other_wind = .true.
               s% kap_rq% use_Type2_opacities = .true.
               s% kap_rq% Zbase = s% initial_z
            end if
            ! if restart during hydro but before CE
            if (s% lxtra(lx_hydro_on)) then
               call star_read_controls(id, 'inlist_hydro_on', ierr)
               s% dt_next = min(1d5, s% dt_next)
               s% dt = min(1d5, s% dt)
               call star_set_u_flag(id, .true., ierr)
            end if
            ! restart during CE
            if (s% lxtra(lx_ce_on)) then
               call star_read_controls(id, 'inlist_ce', ierr)
               r_acc = s% xtra(x_r_acc)  ! do this once to avoid issues calling edep
            end if
            s% use_other_before_struct_burn_mix = .true.
         end if

      end subroutine extras_startup


      real(dp) function envelope_binding_energy(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type(star_info), pointer :: s
         integer :: k
         real(dp) :: envelope_internal_energy
         real(dp) :: envelope_gravitational_energy
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! energy budget integrated from the loaded model
         envelope_internal_energy = 0d0
         do k = 1, s% he_core_k
            envelope_internal_energy = envelope_internal_energy + &
               s% energy(k) * s% dm(k)
         end do
         envelope_gravitational_energy = - dot_product(s% dm_bar(1:s% he_core_k), &
            s% cgrav(1:s% he_core_k) * s% m_grav(1:s% he_core_k) /s% r(1:s% he_core_k))
         envelope_binding_energy = abs(envelope_internal_energy + envelope_gravitational_energy)

      end function envelope_binding_energy


      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         character(len=200) :: fname
         integer :: k, k1, k_peak
         real(dp), pointer :: v(:)
         real(dp) :: tdyn_star, tdyn_core, tdyn
         real(dp) :: m_conv, f_conv
         real(dp) :: vesc, avg_v_div_vesc

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         if (s% u_flag) then
            v => s% u
         else
            v => s% v
         end if

         ! avoid going above MLT velocity limit
         s% max_conv_vel_div_csound = 1d0

         if (abs(log10(abs(s% L_nuc_burn_total * Lsun / s% L(1)))) < 0.005 .and. &
            .not. s% lxtra(lx_found_zams)) then
            if (s% xtra(x_X_central_initial) - s% center_h1 > 0.015d0) then
               write(*,'(/,a,/)') 'Found ZAMS! Turning winds on'
               s% lxtra(lx_found_zams) = .true.
               s% use_other_wind = .true.
               s% kap_rq% use_Type2_opacities = .true.
               s% kap_rq% Zbase = s% initial_z
               s% use_other_before_struct_burn_mix = .true.
               s% job% pgstar_flag = .true.
               ! save a profile, model and photo when reaching ZAMS
               write(fname, '(a22)') 'LOGS/profile_zams.data'
               call star_write_profile_info(id, fname, ierr)
               write(fname, '(a14)') 'model_zams.mod'
               call star_write_model(id, fname, ierr)
               write(fname, '(a17)') 'photos/photo_zams'
               call star_save_for_restart(id, fname, ierr)
            end if
         end if

         
         if (s% lxtra(lx_ce_on)) then
            ! check if during CE we have a shock and reduce max timestep
            k_peak = maxloc(v(1:s% nz), dim=1)
            if (abs(v(k_peak) / s% csound(k_peak)) > 0.9d0 .and. s% m(k_peak) / s% m(1) > 0.99d0) then
               s% max_timestep = 1d3  ! in sec
            end if
            do k=s% nz-1, 1, -1
               vesc = sqrt(2 * s% cgrav(k) * s% m(k) / s% r(k))
               if (v(k) > vesc) then
                  write(*,*) 'found escape velocities'
                  exit
               end if
            end do
            if (k > 1) then
               avg_v_div_vesc = 0d0
               do k1=1, k
                  avg_v_div_vesc = avg_v_div_vesc + &
                     s% dm(k1) * v(k1) / sqrt(2 * s% cgrav(k1) * s% m(k1) / s% r(k1))
               end do
               avg_v_div_vesc = avg_v_div_vesc / (s% m(1) - s% m(k))
               write(*,*) 'average v/vesc', avg_v_div_vesc
            end if
         end if


         ! once the envelope is considered convective, start to model the CE
         ! energy injection, but first turn hydrodynamics on for some time
         ! to have the model settle into a new equilibrium condition
         if (s% lxtra(lx_hydro_on) .and. .not. s% lxtra(lx_ce_on)) then
            write(*,'(/,a24,2x,1pd8.2,2x,a3)') 'time left to turn CE on:', &
               abs(s% time - (s% xtra(x_time_start_ce) + num_dyn * s% xtra(x_dyn_time))) / secyer, 'yrs'
            if (s% time > s% xtra(x_time_start_ce) + num_dyn * s% xtra(x_dyn_time) .and. &
               .not. s% lxtra(lx_ce_on)) then
               write(*,'(/,a,/)') 'model relaxed after hydro on. turning CE on'
               s% lxtra(lx_ce_on) = .true.
               call star_read_controls(id, 'inlist_ce', ierr)
               call star_set_conv_vel_flag(id, .false., ierr)
               s% use_other_before_struct_burn_mix = .true.
               s% xtra(x_r_acc) = fraction_of_r_donor * s% r(1)
               r_acc = s% xtra(x_r_acc)  ! do this once to avoid issues calling edep
            end if
         end if


         ! search for convective envelope condition
         if (.not. s% lxtra(lx_have_convective_envelope)) then
            if (s% center_h1 < 1d-4) then
               m_conv = 0d0
               do k=1, s% he_core_k
                  if (s% mixing_type(k) == convective_mixing) m_conv = m_conv + s% dm(k)
               end do
               m_conv = m_conv / Msun
               f_conv = m_conv / s% star_mass
               write(*,'(a40, f26.16)') 'convective env mass', m_conv
               write(*,'(a40, f26.16)') 'fract of conv env', f_conv
               if (f_conv > 0.2d0) then
                  write(*,'(/,a,/)') 'convection dominates in more than 20% of the mass of entire star. turning hydro on'
                  s% lxtra(lx_have_convective_envelope) = .true.
                  s% lxtra(lx_hydro_on) = .true.
                  call star_read_controls(id, 'inlist_hydro_on', ierr)
                  s% dt_next = min(1d2, s% dt_next)
                  s% dt = min(1d2, s% dt)
                  call star_set_u_flag(id, .true., ierr)
                  ! force initial velocities to zero to prevent issues in outer layers
                  s% xh(s% i_u,:) = 0d0
                  s% xh_old(s% i_u,:) = 0d0
                  s% generations = 1
                  ! eval dynamic timescale of envelope
                  tdyn_star = 1d0 / sqrt(s% m(1) / (4d0/3d0 * pi * pow3(s% r(1))) * standard_cgrav)
                  tdyn_core = 1d0 / sqrt(s% m(s% he_core_k) / (4d0/3d0 * pi * pow3(s% r(s% he_core_k))) * standard_cgrav)
                  tdyn = tdyn_star - tdyn_core  ! in sec
                  s% xtra(x_dyn_time) = tdyn
                  s% xtra(x_time_start_ce) = s% time  ! in sec
                  write(*,'(a14,2x,1pd8.2,2x,a3)') 'will be on for', num_dyn * s% xtra(x_dyn_time) / secyer, 'yrs'
                  ! get binding energy of envelope
                  s% xtra(x_initial_energy) = envelope_binding_energy(id, ierr)
                  if (ierr /= 0) return
                  if (dbg) write(*,'(/,a40, 1pd26.16,/)') 'energy_budget', s% xtra(x_initial_energy)
                  ! save a profile, model and photo when reaching convective envelope condition
                  write(fname, '(a37)') 'LOGS/profile_convective_envelope.data'
                  call star_write_profile_info(id, fname, ierr)
                  write(fname, '(a29)') 'model_convective_envelope.mod'
                  call star_write_model(id, fname, ierr)
                  write(fname, '(a21)') 'photos/photo_hydro_on'
                  call star_save_for_restart(id, fname, ierr)
               end if
            end if
         end if

         if (s% use_other_before_struct_burn_mix) then
            call before_struct_burn_mix_for_edep(id, s% dt, extras_start_step)
         end if

         extras_start_step = 0

         if (s% lxtra(lx_hydro_on)) then
            !!write(*,*) 'use_ODE_var_eqn_pairing', s% use_ODE_var_eqn_pairing
            write(*,*) 'use_dPrad', s% use_dPrad_dm_form_of_T_gradient_eqn
            write(*,*) 'use_dedt', s% use_dedt_form_of_energy_eqn
            write(*,*) 'use_moomentum', s% use_momentum_outer_BC
            write(*,*) 'use_compression', s% use_compression_outer_BC
         end if

   
      end function extras_start_step


      subroutine before_struct_burn_mix_for_edep(id, dt, res)
         integer, intent(in) :: id
         real(dp), intent(in) :: dt
         integer, intent(out) :: res ! keep_going, redo, retry, backup, terminate
         integer :: ierr
         type (star_info), pointer :: s
         real(dp), pointer :: v(:)
         integer :: k_peak

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         if (s% u_flag) then
            v => s% u
         else
            v => s% v
         end if

         if (s% lxtra(lx_hydro_on)) then
            call star_read_controls(id, 'inlist_hydro_on', ierr)
            if (s% lxtra(lx_ce_on)) then
               call star_read_controls(id, 'inlist_ce', ierr)
               k_peak = maxloc(v(1:s% nz), dim=1)
               if (abs(v(k_peak) / s% csound(k_peak)) > 0.9d0 .and. s% m(k_peak) / s% m(1) > 0.99d0) then
                  write(*,*) 'shock near surface. changing BC'
                  s% use_momentum_outer_BC = .true.
                  s% use_compression_outer_BC = .false.
               else
                  s% use_compression_outer_BC = .true.
                  s% use_momentum_outer_BC = .false.
               end if
            end if
         else
            call star_read_controls(id, 'inlist_hydro_off', ierr)
            s% kap_rq% Zbase = s% initial_z
         end if

         ! winds off during hydro
         if (s% lxtra(lx_hydro_on)) then
            s% use_other_wind = .false.
            s% tau_factor = 1
            s% force_tau_factor = 1
            s% tau_base = 1d-4
         end if

         ! according to ppisn test_case this flag can be override when reading inlists
         s% use_other_before_struct_burn_mix = .true.
         
         res = keep_going

      end subroutine before_struct_burn_mix_for_edep


      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         extras_check_model = keep_going

      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         how_many_extra_history_columns = 14

      end function how_many_extra_history_columns

      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k_peak
         real(dp) :: vesc_at_peak
         real(dp), pointer :: v(:)

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! s% xtra(x_r_acc) = r_acc
         ! s% xtra(x_dr) = dr
         ! s% xtra(x_bondi_radius) = bondi_radius
         ! s% xtra(x_de) = de
         ! s% xtra(x_v_circ) = v_circ
         ! s% xtra(x_cs) = cs
         ! s% xtra(x_penetration_coeff) = penetration_coeff
         ! s% xtra(x_eps_rho) = eps_rho
         ! s% xtra(x_v_esc) = v_esc
         ! s% xtra(x_r_at_v_esc) = r_at_v_esc
         ! s% xtra(x_v_at_v_esc) = v_at_v_esc

         names(1) = 'r_acc'
         names(2) = 'dr'
         names(3) = 'bondi_radius'
         names(4) = 'de'
         names(5) = 'orbital_velocity'
         names(6) = 'cs'
         names(7) = 'mach_number'
         names(8) = 'eps_rho'
         names(9) = 'v_esc'
         names(10) = 'r_at_v_esc'
         names(11) = 'v_at_v_esc'
         vals(1) = s% xtra(x_r_acc) / Rsun
         vals(2) = s% xtra(x_dr) / Rsun
         vals(3) = s% xtra(x_bondi_radius) / Rsun
         vals(4) = s% xtra(x_de)
         vals(5) = s% xtra(x_v_circ) / 1d5
         vals(6) = s% xtra(x_cs)
         vals(7) = s% xtra(x_v_circ) / s% xtra(x_cs)
         vals(8) = s% xtra(x_eps_rho)
         vals(9) = s% xtra(x_v_esc)
         vals(10) = s% xtra(x_r_at_v_esc)
         vals(11) = s% xtra(x_v_at_v_esc)

         names(12) = 'r_at_peak'
         names(13) = 'v_div_cs_at_peak'
         names(14) = 'vesc_div_cs_at_peak'
         if (s% lxtra(lx_hydro_on)) then
            if (s% u_flag) then
               v => s% u
            else
               v => s% v
            end if
            k_peak = maxloc(v(1:s% nz), dim=1)
            vesc_at_peak = sqrt(2 * s% cgrav(k_peak) * s% m(k_peak) / s% r(k_peak))
            vals(12) = s% r(k_peak) / Rsun
            vals(13) = v(k_peak) / s% csound(k_peak)
            vals(14) = vesc_at_peak / s% csound(k_peak)
         else
            vals(12) = 0d0
            vals(13) = 0d0
            vals(14) = 0d0
         end if

      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         how_many_extra_profile_columns = 6

      end function how_many_extra_profile_columns

      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp), pointer :: v(:)
         integer :: k
         real(dp) :: v_esc

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         if (s% u_flag) then
            v => s% u
         else
            v => s% v
         end if

         names(1) = 'v_div_vesc'
         names(2) = 'specific_e_grav'
         names(3) = 'specific_e_kin'
         names(4) = 'specific_e_th'
         names(5) = 'total_specific_e'
         names(6) = 'mlt_vc'
         do k=1, nz
            v_esc = sqrt(2d0 * s% cgrav(k) * s% m(k) / s% r(k))
            vals(k,1) = v(k) / v_esc
            vals(k,2) = -s% cgrav(k) * s% m(k) / s% r(k)
            vals(k,3) = 0.5d0 * v(k) * v(k)
            vals(k,4) = two_thirds * avo * kerg* s% T(k) / (2 * s% mu(k) * s% rho(k)) &
               + crad * pow4(s% T(k)) / s% rho(k)
            vals(k,5) = vals(k,2) + vals(k,3) + vals(k,4)
            vals(k,6) = s% mlt_vc(k)
         end do

      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         how_many_extra_history_header_items = 0

      end function how_many_extra_history_header_items

      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr

         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         how_many_extra_profile_header_items = 0

      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr

         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         character (len=200) :: fname
         integer :: k, nz
         logical :: going_to_remove_surface
         real(dp) :: q, rlobe
         real(dp) :: v_esc
         real(dp), pointer :: v(:)

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         extras_finish_step = keep_going


         ! if not using energy injection, do not evaluate anything
         if (.not. s% use_other_energy) return
         
         nz = s% nz
         going_to_remove_surface = .false.
         if (s% u_flag) then
            v => s% u
         else
            v => s% v
         end if

         ! warn if we have reached the initial energy budget
         if (s% cumulative_extra_heating >= s% xtra(x_initial_energy)) &
            write(*,'(/,a,/)') 'WARNING: have injected more than initial binding energy'

         ! check for the overflow of the he-core, in that case stop
         !q = s% m(s% he_core_k) / m_acc
         !rlobe = r_acc * 0.49d0 * q * q / (0.6d0 * q * q + log1p(q))
         !if (s% he_core_radius * Rsun >= rlobe) then
         !   write(*,'(a)') 'terminate due to he-core RLOF'
         !   write(*,'(a40, f26.16)') 'core radius', s% he_core_radius
         !   write(*,'(a40, f26.16)') 'rlobe', rlobe / Rsun
         !   write(*,'(a)')
         !   extras_finish_step = terminate
         !   return
         !end if

         ! update value of r_acc
         r_acc = max(s% xtra(x_r_acc) - s% xtra(x_dr), 0d0)
         write(*,'(a40, f26.16)') 'r_acc updated', r_acc / Rsun
         if (r_acc > s% r(1)) then
            write(*,'(a)') 'forcing accretor to be inside donor'
            r_acc = max(s% r(1) - s% xtra(x_dr), 0d0)
         end if


         ! check for escape velocities near surface 
         do k=nz, 1, -1
            v_esc = sqrt(2d0 * s% cgrav(k) * s% m(k) / s% r(k))
            if (v(k) > v_esc .and. s% q(k) > 0.9d0) then
               write(*,'(a)') 'found place where velocity is bigger than local escape velocity'
               write(*,'(a40, f26.16)') 'v_esc', v_esc
               write(*,'(a40, f26.16)') 'R', s% r(k) / Rsun
               if (remove_surface) going_to_remove_surface = .true.
               s% xtra(x_v_esc) = v_esc
               s% xtra(x_r_at_v_esc) = s% r(k) / Rsun
               s% xtra(x_v_at_v_esc) = v(k)
               exit
            end if
         end do

         if (going_to_remove_surface) then
            write(*,'(a40, i14, /)') 'removing surface above cell k', k
            call star_remove_surface_at_cell_k(id, k, ierr)
            if (ierr /= 0) then
               stop 'problem removing surface'
               return
            end if
         end if

         ! also check for low densities and remove them
         do k=nz, 1, -1
            if (s% r(k) / Rsun > 1d4) then
               write(*,'(a40, i14, /)') 'found place where radius is above 10**4', k
               call star_remove_surface_at_cell_k(id, k, ierr)
               exit
            end if
         end do

      end function extras_finish_step


      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call test_suite_after_evolve(s, ierr)

      end subroutine extras_after_evolve


      subroutine alloc_extra_info(s)
         integer, parameter :: extra_info_alloc = 1
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_alloc)
      end subroutine alloc_extra_info


      subroutine unpack_extra_info(s)
         integer, parameter :: extra_info_get = 2
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_get)
      end subroutine unpack_extra_info


      subroutine store_extra_info(s)
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_put)
      end subroutine store_extra_info


      subroutine move_extra_info(s,op)
         integer, parameter :: extra_info_alloc = 1
         integer, parameter :: extra_info_get = 2
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         integer, intent(in) :: op

         integer :: i, j, num_ints, num_dbls, ierr

         i = 0
         ! call move_int or move_flg
         call move_int(redo_count)
         num_ints = i

         i = 0
         num_dbls = i

         if (op /= extra_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return

         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
         end if

         contains

         subroutine move_dbl(dbl)
            double precision :: dbl
            i = i+1
            select case (op)
            case (extra_info_get)
               dbl = s% extra_work(i)
            case (extra_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl

         subroutine move_int(int)
            integer :: int
            i = i+1
            select case (op)
            case (extra_info_get)
               int = s% extra_iwork(i)
            case (extra_info_put)
               s% extra_iwork(i) = int
            end select
         end subroutine move_int

         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (extra_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (extra_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg

      end subroutine move_extra_info

      end module run_star_extras
      
