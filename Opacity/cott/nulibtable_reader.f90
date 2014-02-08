!-*-f90-*-
subroutine nulibtable_reader(filenames)
  
  use nulibtable
  use hdf5
  implicit none

  character(len=512), dimension(3) :: filenames

  integer :: inu
  character(len=512) :: filename
  
  !H5 stuff
  integer :: error,rank,cerror
  integer(HID_T) :: file_id,dset_id,dspace_id
  integer(HSIZE_T) :: dims1(1), dims2(2), dims3(3), dims4(4), dims5(5)!, etc....

  !local
  real*8, allocatable :: nulibtable_temp(:,:,:,:,:)
  real*8 :: timestamp
  integer :: startindex,endindex
  integer :: i,ngrp, nrho, nye, ntemp
  type(neutrino), pointer :: neutp

  cerror = 0

  nulibtable_number_species = 0

  do inu=1,3

     if (inu .eq. 1) then
        neutp => electron_neut
     else if (inu .eq. 2) then
        neutp => electron_antineut
     else
        neutp => muon_neut
     end if

     filename = filenames(inu)

     !open HDF5 file, given filename                                                                                             
     call h5open_f(error)
     if (error.ne.0) then
        stop "Error reading in nulib file"
     endif
     
     call h5fopen_f(trim(adjustl(filename)), &
          H5F_ACC_RDONLY_F,file_id,error)
     if (error.ne.0) then
        write(*,*) trim(adjustl(filename))
        stop "Error reading in nulib table"
     endif
     
     !write scalars (rank=1, dims1(1) = 1)
     rank = 1
     dims1(1) = 1
     
     !first lets read number of species and number of groups and number of rho,temp, and ye points, also timestamp
     nulibtable_number_species = nulibtable_number_species + 1

     call h5dopen_f(file_id, "timestamp",dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, timestamp, dims1, error)
     call h5dclose_f(dset_id, error)
     cerror = cerror + error  

     call h5dopen_f(file_id, "nrho",dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nrho, dims1, error)
     call h5dclose_f(dset_id, error)
     cerror = cerror + error

     if (inu .eq. 1) then
        nulibtable_nrho = nrho
     else
        if (nrho .ne. nulibtable_nrho) then
           print *, "Wrong nrho: ", inu, nrho, nulibtable_nrho
           stop
        end if
     end if
     
     call h5dopen_f(file_id, "ntemp",dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER, ntemp, dims1, error)
     call h5dclose_f(dset_id, error)
     cerror = cerror + error

     if (inu .eq. 1) then
        nulibtable_ntemp = ntemp
     else
        if (ntemp .ne. nulibtable_ntemp) then
           print *, "Wrong ntemp: ", inu, ntemp, nulibtable_ntemp
           stop
        end if
     end if

     call h5dopen_f(file_id, "nye",dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nye, dims1, error)
     call h5dclose_f(dset_id, error)
     cerror = cerror + error

     if (inu .eq. 1) then
        nulibtable_nye = nye
     else
        if (nye .ne. nulibtable_nye) then
           print *, "Wrong nye: ", inu, nye, nulibtable_nye
           stop
        end if
     end if

     call h5dopen_f(file_id, "number_groups",dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER, ngrp, dims1, error)
     call h5dclose_f(dset_id, error)
     cerror = cerror + error
     
     call init_neutrino(neutp, ngrp)

     !lets also read neutrino energies, bin bottoms,tops and widths
     rank = 1
     dims1(1) = ngrp
     call h5dopen_f(file_id, "neutrino_energies",dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,neutp%energies, dims1, error)
     call h5dclose_f(dset_id, error)
     cerror = cerror + error      
  
     !calculate inverse energies
     neutp%inv_energies = 1.0d0/neutp%energies

     call h5dopen_f(file_id, "bin_widths", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,neutp%ewidths, dims1, error)
     call h5dclose_f(dset_id, error)
     cerror = cerror + error   

     call h5dopen_f(file_id, "bin_bottom", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,neutp%ebottom, dims1, error)
     call h5dclose_f(dset_id, error)
     cerror = cerror + error   

     call h5dopen_f(file_id, "bin_top", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,neutp%etop, dims1, error)
     call h5dclose_f(dset_id, error)
     cerror = cerror + error   

     if (inu .eq. 1) then
        allocate(nulibtable_logrho(nulibtable_nrho))
        allocate(nulibtable_logtemp(nulibtable_ntemp))
        allocate(nulibtable_ye(nulibtable_nye))
    
        rank = 1
        dims1(1) = nulibtable_nrho
        call h5dopen_f(file_id, "rho_points", dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,nulibtable_logrho, dims1, error)
        call h5dclose_f(dset_id, error)
        cerror = cerror + error     

        rank = 1
        dims1(1) = nulibtable_ntemp
        call h5dopen_f(file_id, "temp_points", dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,nulibtable_logtemp, dims1, error)
        call h5dclose_f(dset_id, error)
        cerror = cerror + error  
        
        rank = 1
        dims1(1) = nulibtable_nye
        call h5dopen_f(file_id, "ye_points", dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,nulibtable_ye, dims1, error)
        call h5dclose_f(dset_id, error)
        cerror = cerror + error   

        nulibtable_logrho = log10(nulibtable_logrho)
        nulibtable_logrho_min = nulibtable_logrho(1)
        nulibtable_logrho_max = nulibtable_logrho(nulibtable_nrho)
        nulibtable_logtemp = log10(nulibtable_logtemp) 
        nulibtable_logtemp_min = nulibtable_logtemp(1)
        nulibtable_logtemp_max = nulibtable_logtemp(nulibtable_ntemp)
        nulibtable_ye_min = nulibtable_ye(1)
        nulibtable_ye_max = nulibtable_ye(nulibtable_nye)

     end if
  
     !now read three tables
     rank = 5
     dims5(1) = nulibtable_nrho
     dims5(2) = nulibtable_ntemp
     dims5(3) = nulibtable_nye
     dims5(4) = 1
     dims5(5) = ngrp  

     allocate(nulibtable_temp(nulibtable_nrho,nulibtable_ntemp, &
          nulibtable_nye,1,ngrp))  

     call h5dopen_f(file_id, "emissivities", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,nulibtable_temp, dims5, error)
     call h5dclose_f(dset_id, error)
     cerror = cerror + error   

     do i=1,ngrp
        neutp%emissivities(:,:,:,i) = log10(max(1.0d-30,&
             nulibtable_temp(:,:,:,1,i)*neutp%ewidths(i)))
     end do

     call h5dopen_f(file_id, "absorption_opacity",dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,nulibtable_temp, dims5, error)
     call h5dclose_f(dset_id, error)
     cerror = cerror + error   

     neutp%absopacity = log10(nulibtable_temp(:,:,:,1,:))

     call h5dopen_f(file_id, "scattering_opacity", dset_id, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,nulibtable_temp, dims5, error)
     call h5dclose_f(dset_id, error)
     cerror = cerror + error   

     neutp%scatopacity = log10(nulibtable_temp(:,:,:,1,:))

     deallocate(nulibtable_temp)

     !must close h5 files, check for error
     if (cerror.ne.0) then
        write (*,*) "We have errors on reading HDF5 file", cerror
        stop
     endif
     
     call h5fclose_f(file_id,error)
     call h5close_f(error)
     
     nullify(neutp)
     
  end do

end subroutine nulibtable_reader
