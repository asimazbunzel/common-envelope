      ! orbital velocity (circular)
      real(dp) function orbital_velocity(m1, r)
         real(dp), intent(in) :: m1, r

         orbital_velocity = sqrt(standard_cgrav * m1 / r)

      end function orbital_velocity


      ! bondi accretion radius
      real(dp) function bondi_accretion_radius(m, cs, v)
         real(dp), intent(in) :: m, cs, v

         if (v < cs) then
            bondi_accretion_radius = 2d0 * standard_cgrav * m / (pow2(v) + pow2(cs))
         else
            bondi_accretion_radius = 2d0 * standard_cgrav * m / pow2(v)
         end if

      end function bondi_accretion_radius


      ! get cell values matching (r_acc - bondi_radius, r_acc + bondi_radius)
      subroutine get_influence_cell_k(nz, r_star, r_acc, bondi_radius, k_bot, k_top)
         integer, intent(in) :: nz
         real(dp), intent(in) :: r_star(nz)
         real(dp), intent(in) :: r_acc, bondi_radius
         integer, intent(out) :: k_bot, k_top
         integer :: k

         k_bot = 2
         do k = 2, nz
            if (r_star(k) < r_acc - bondi_radius) then
               k_bot = k
               exit
            end if
         end do

         k_top = 1
         do k = 1, nz
            if (r_star(k) < r_acc + bondi_radius) then
               k_top = k
               exit
            end if
         end do

      end subroutine get_influence_cell_k


      ! penetration depth function
      real(dp) function eval_penetration(bondi_radius, Rstar, r_wd) result(penetration_coeff)
         real(dp), intent(in) :: bondi_radius, Rstar, r_wd

         penetration_coeff = bondi_radius + (Rstar - r_wd)

         if (penetration_coeff < 0d0) stop 'penetration_coeff < 0'

      end function eval_penetration


      ! engulfed area of r_wd with an accretion radius from the HLA theory
      real(dp) function engulfed_area(penetration_coeff, bondi_radius) result(area)
         real(dp), intent(in) :: penetration_coeff, bondi_radius
         real(dp) :: y, alpha

         if (penetration_coeff < bondi_radius) then
            y = bondi_radius - penetration_coeff
            alpha = acos(y / bondi_radius)
            area = bondi_radius * (bondi_radius * alpha - y * sin(alpha))
         else
            y = penetration_coeff - bondi_radius
            alpha = acos(y / bondi_radius)
            area = pi * pow2(bondi_radius) - bondi_radius * (bondi_radius * alpha - y * sin(alpha))
         end if

      end function engulfed_area


      ! gravitational drag force on the wd
      ! for it, we use HLA theory which, after equating dot_Eorb from Kepler
      ! to dot_Eorb = F_HLA * v_circ gives the new position of the WD:
      ! dr = (2 * area * rho * sqrt(G * m_acc * r_acc) / mr) * dt
      subroutine eval_drag_force(mr, m_acc, area, rho, dt, r, de, dr)
         real(dp), intent(in) :: mr, m_acc, area, rho, dt, r
         real(dp), intent(out) :: de, dr
         real(dp) :: cdr, cde, csi

         ! csi factor to increase in supersonic motion (Shima et al., 1985)
         csi = 4d0

         cdr = csi * 2d0 * area * rho * sqrt(standard_cgrav * m_acc) / mr
         dr = cdr * sqrt(r) * dt

         cde = csi * 0.5d0 * standard_cgrav * m_acc * mr / pow2(r)
         de = cde * dr

      end subroutine eval_drag_force


      ! correction to the classical gravitational drag force in a supersonic regime
      subroutine macleod_fit(s, k, accretion_radius, eps_rho, fit)
         integer, intent(in) :: k
         real(dp), intent(in) :: accretion_radius
         real(dp), intent(out) :: eps_rho
         real(dp), intent(out) :: fit
         type(star_info), pointer :: s
         real(dp) :: dr_drho, density_scale_height

         dr_drho = (s% r(k+1) - s% r(k)) / (s% rho_face(k+1) - s% rho_face(k))
         density_scale_height = - s% rho_face(k) * dr_drho

         eps_rho = accretion_radius / density_scale_height

         fit = 1.91791946d0 - 1.52814698d0 * eps_rho + 0.75992092d0 * pow2(eps_rho)

         if (fit < 0d0) stop 'f_m below 0'

      end subroutine macleod_fit
