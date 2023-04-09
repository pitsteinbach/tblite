! Copyright (C) 2021 Christoph Bannwarth
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

!> compute atomic response function and effective H-L gap
! nat:  # of atoms
! nao:  # of AOs/MOs (assuming spherical basis here)
! focc: (fractional) occupation numbers
! eps:  orbital energies
! aoat: ao-to-orbital indexer
! C: MO coefficient
! S: overlap matrix
! response: atomic response (measures level availability on atom)
! egap: effective atomic H-L gap
! chempot: efffective Fermi level per atom (as 0.5*(eh+el)
! ehoao: highest occ. atomic orbital
! eluao: lowest virt. atomic orbital
subroutine atomic_frontier_orbitals(nat,nao,focc,eps,aoat,C,S,response,egap,chempot,ehoao,eluao)
   use mctc_env, only : wp,sp
   !use xtb_mctc_accuracy, only : wp, sp
   use mctc_io_convert, only : autoev
   

   integer, intent(in)    :: nat,nao,aoat(nao)
   real(wp),intent(in)    :: focc(nao),eps(nao)
   real(wp),intent(in)    :: C(nao,nao),S(nao,nao)

   real(wp), intent(out)  :: ehoao(nat),eluao(nat),egap(nat),response(nat),chempot(nat)
   integer  :: i,j,k,jj,kk,m
   real(wp),allocatable :: po(:,:),pv(:,:)
   real(wp) :: occ,tmp,tmp2,ps,virt,tmp3,weight,osite,vsite
   real(wp), parameter :: damp=0.5_wp ! damping in response function (in eV)
   allocate(po(nat,nao),pv(nat,nao))

   po=0.0_wp
   pv=0.0_wp
   
   
   ! collect the atomic MO populations
   ! we collect occ & virt separately (and allow for fractional occupation)
   do i=1,nao
      ! occ part
      occ=0.50_wp*focc(i)
      if(occ.gt.0.0001_wp) then 
        do j=1,nao
            jj=aoat(j)
            do k=1,j-1
               kk=aoat(k)
               ps=s(k,j)*C(j,i)*C(k,i)*occ
               po(kk,i)=po(kk,i)+ps
               po(jj,i)=po(jj,i)+ps
            enddo
            ps=s(j,j)*C(j,i)*C(k,i)*occ
            po(jj,i)=po(jj,i)+ps
        enddo
      endif
   enddo

   do i=1,nao
      !virt part
      virt = 0.50_wp *(2.0_wp - focc(i))
      if(virt.gt.0.0001_wp) then 
        do j=1,nao
            jj=aoat(j)
            do k=1,j-1
               kk=aoat(k)
               ps=s(k,j)*C(j,i)*C(k,i)*virt
               pv(kk,i)=pv(kk,i)+ps
               pv(jj,i)=pv(jj,i)+ps
            enddo
            ps=s(j,j)*C(j,i)*C(k,i)*virt
            pv(jj,i)=pv(jj,i)+ps
        enddo
      endif
          
   enddo
!   write(*,*) (sum(po(i,:)),i=1,nat) ! DEBUG
!   write(*,*) (sum(pv(i,:)),i=1,nat) ! DEBUG

   ! now accumulate the atomic H-L gaps and Fermi levels (we approximate it as 0.5*(eH + eL)) 
   ! we make use of an atomic response-type weightinhg: 
   ! chempot=\sum_ai * wA_ai * 0.5 * (e_a + e_i) 
   ! here wA_ai is the variable "response" computed as: wA_ai=  [p_i p_a / ( (e_a - e_i)**2 + damp**2) ]/ [\sum_ia p_j p_b / ((e_b - e_j)**2 + damp**2)]
   ! the "p"s are the MO densities, for gaps we accumulate the regularized inverse and later on invert it again
   response=0.0_wp
   chempot=0.0_wp
   egap=0.0_wp
   ! go through occ (including fractionally occupied)
   do i=1,nao
      occ =0.50_wp * focc(i)
      if(occ.gt.0.0001_wp) then
        osite=sum(po(:,i))/maxval(po(:,i)) ! not used, but roughly gives the number of centers for the canonical MOs
        !virt part  (including fractionally occupied)
        do j=1,nao 
           virt = 0.5_wp *(2.0_wp - focc(j))
           if(virt.gt.0.0001_wp) then
              vsite=sum(pv(:,j))/maxval(pv(:,j)) ! not used, but roughly gives the number of centers for the canonical MOs
              tmp=1.0_wp/((eps(j)-eps(i))**2 + damp**2)
              tmp2=0.5_wp * (eps(j)+eps(i))
              tmp3=1.0_wp /(eps(j)-eps(i) +damp )
              ! compute atomic response 
              do m=1,nat
                  weight= po(m,i)*pv(m,j)*tmp
                  response(m)=response(m)+weight
                  chempot(m)=chempot(m)+tmp2*weight
                  egap(m)=egap(m)+weight*tmp3
              enddo ! at
          endif ! if virt 
        enddo !  virt  
      endif ! if occ
          
   enddo ! occ

   ehoao=0.0_wp
   eluao=0.0_wp

   write(*,*)
   write(*,*) "  -------------------------"
   write(*,*) "  atomic frontier MO info "
   write(*,*) "  -------------------------"
   write(*,*) "  atom   response (a.u.)   gap (eV)  chem.pot (eV)  HOAO (eV)    LUAO (eV)" 
   ! now get the atomic frontier orbitals
   tmp=0.0_wp
   
   do m=1,nat
      egap(m)=egap(m)/response(m)
      egap(m)=1.0_wp/egap(m)-damp
      chempot(m)=chempot(m)/response(m)
      ehoao(m)=(chempot(m)-0.5_wp*egap(m))
      eluao(m)=(chempot(m)+0.5_wp*egap(m))
      write(*,'(3x,i6,5x,f12.4,5x,f8.4,5x,f8.4,5x,f8.4,5x,f8.4)') m,response(m)*(autoev**2),egap(m),chempot(m),ehoao(m),eluao(m)

   enddo
   write(*,*)
   deallocate(po,pv)
end subroutine atomic_frontier_orbitals

