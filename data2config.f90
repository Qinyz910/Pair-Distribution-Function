!--------------------------------------------------------------------------
!
! data2config
!
! Generates atomic configuration for large-scale simulations starting from
! crystal structure data or some particular configuration file formats.
!
! Mostly designed to support RMCprofile, but is also useful for DL_POLY
!
! Written by Martin Dove.
! DLPOLY part co-written by Mike Andrew
! Debugging and testing by Matt Tucker
!
! See module version_data for the version number and date
!
!--------------------------------------------------------------------------
!
! Notes on how to add additional functionality
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! If I want to add support for another input file format ...
! 1. Read file type from filename extension in the subroutine get_arguments.
! 2. Set up file type in character variable cinput.
! 3. Ensure that if I can read from the same type of file that the logical
!    flag is set as .false..
!
! If I want to add support for another output file format ...
! 1. Create a logical flag for this, declare it in the module arguments, and
!    set it in the subroutine get_arguments.
! 2. Ensure that I am not reading from the same file type, by setting the argument
!    for this flag as .false. in the subroutine get_arguments.
! 3. Create the input and output subroutines.
!
!--------------------------------------------------------------------------

module version_data

    integer, parameter :: nversions = 29
    character(len=120) :: version_history(nversions)
    character(len=5) :: version_number(nversions)
    character(len=20) :: version_date(nversions)
    
    data version_number(1) / '3.10' /
    data version_number(2) / '3.11' /
    data version_number(3) / '3.12' /
    data version_number(4) / '3.13' /
    data version_number(5) / '3.14' /
    data version_number(6) / '3.15' /
    data version_number(7) / '3.16' /
    data version_number(8) / '3.17' /
    data version_number(9) / '3.18' /
    data version_number(10) / '3.19' /
    data version_number(11) / '4.0' /
    data version_number(12) / '4.1' /
    data version_number(13) / '4.2' /
    data version_number(14) / '4.3' /
    data version_number(15) / '4.4' /
    data version_number(16) / '4.5' /
    data version_number(17) / '4.5b' /
    data version_number(18) / '4.5c' /
    data version_number(19) / '4.6a' /
    data version_number(20) / '4.6b' /
    data version_number(21) / '4.7' /
    data version_number(22) / '4.8' /
    data version_number(23) / '4.9' /
    data version_number(24) / '5.0' /
    data version_number(25) / '5.1' /
    data version_number(26) / '5.2' /
    data version_number(27) / '5.3' /
    data version_number(28) / '5.4' /
    data version_number(29) / '5.5' /
    
    data version_date(1) / 'May 2016' /
    data version_date(2) / 'June 2016' /
    data version_date(3) / 'July 2016' /
    data version_date(4) / 'July 2016' /
    data version_date(5) / 'August 2016' /
    data version_date(6) / 'September 2016' /
    data version_date(7) / 'December 2016' /
    data version_date(8) / 'January 1 2017' /
    data version_date(9) / 'January 19 2017' /
    data version_date(10) / 'February 15 2017' /
    data version_date(11) / 'April 23 2017' /
    data version_date(12) / 'June 1 2017' /
    data version_date(13) / 'June 11 2017' /
    data version_date(14) / 'July 29 2017' /
    data version_date(15) / 'August 21 2017' /
    data version_date(16) / 'September 9 2017' /
    data version_date(17) / 'October 21 2017' /
    data version_date(18) / 'October 28 2017' /
    data version_date(19) / 'November 29 2017' /
    data version_date(20) / 'December 26 2017' /
    data version_date(21) / 'May 3 2018' /
    data version_date(22) / 'May 6 2018' /
    data version_date(23) / 'August 1 2018' /
    data version_date(24) / 'December 23 2018' /
    data version_date(25) / 'February 19 2019' /
    data version_date(26) / 'May 3 2019' /
    data version_date(27) / 'Janary 23 2020' /
    data version_date(28) / 'September 9 2020' /
    data version_date(29) / 'December 7 2020' /

    data version_history(1) / 'Create gasp configuration' /
    data version_history(2) / 'Read from atomeye .cfg file' /
    data version_history(3) / 'Records version history with -history command' /
    data version_history(4) / 'Prints cumulative n(r) function from the PDF' /
    data version_history(6) / 'Takes a list (AtomEye brief format) to define origins for analysis' /
    data version_history(7) / 'Allows outputs to be selected as fractional or orthogonal coordinates' /
    data version_history(8) / 'Converts DL_POLY HISTORY file to individual files' /
    data version_history(9) / 'Added construction of DL_POLY PDF from RDFDAT file' /
    data version_history(10) / 'Added XTL file output' /
    data version_history(11) / 'Added .cell file, element relabelling, occupancies' /
    data version_history(12) / 'Added write castep cell file' /
    data version_history(13) / 'Added CrystalMaker-driven bondlength writing' /
    data version_history(14) / 'Added random molecule orientations' /
    data version_history(15) / 'Added read from SHELX/CASTEP .res file' /
    data version_history(16) / 'Added dipole analysis' /
    data version_history(17) / 'Rebuilt random rotations from v4.3' /
    data version_history(18) / 'Fixed a problem with creating vacancies' /
    data version_history(19) / 'Applying -limits data to some configuration writing' /
    data version_history(20) / 'Fixed issues of random orientaions and periodic boundaries' /
    data version_history(21) / '-compare option to compare two rmc6f configuration files' /
    data version_history(22) / '-align option to align the structure with a supplied reference rmc6f configuration file' /
    data version_history(23) / '-oneunitcell option to perform analysis using one unit cell only' /
    data version_history(24) / '-nanoscale option to scale a nanoparticle with molecules' /
    data version_history(25) / 'Input www file and use -silica option to generate silica' /
    data version_history(26) / 'Calculate orientational variables and Kubic harmonic functions' /
    data version_history(27) / 'Added errors to averages, analysis of DLPOLY configurations from FIELD file' /
    data version_history(28) / 'Output www-format' /
    data version_history(29) / '-rmcanal, generate analysis using RMC molecules files' /

end module version_data

!===============================================================================

module channel_numbers

    integer,parameter:: ic = 12
    integer,parameter:: main_output = 13
    integer,parameter:: icssr = 14
    integer,parameter:: idlpoly = 15
    integer,parameter:: irmc3 = 16
    integer,parameter:: irmc6f = 17
    integer,parameter:: irmc6o = 18
    integer,parameter:: ihis6f = 19
    integer,parameter:: ihis6o = 20
    integer,parameter:: icrystal = 21
    integer,parameter:: icrush = 22
    integer,parameter:: imag = 23
    integer,parameter:: ifield = 24
    integer,parameter:: ifieldout = 25
    integer,parameter:: ianal = 26
    integer,parameter:: ibondfile = 27
    integer,parameter:: iexp = 28
    integer,parameter:: ilst = 29
    integer,parameter:: iinst = 30
    integer,parameter:: iback = 31
    integer,parameter:: ibragg = 32
    integer,parameter:: icif = 33
    integer,parameter:: ipdfn = 34
    integer,parameter:: ipdfx = 35
    integer,parameter:: igulp = 36
    integer,parameter:: iatomeye = 37
    integer,parameter:: idw = 38
    integer,parameter:: idwout = 39
    integer,parameter:: itemplate = 40
    integer,parameter:: icmtx = 41
    integer,parameter:: imcfg = 42
    integer,parameter:: ireserve = 43  ! reserved for use with icmtx
    integer,parameter:: igasp = 44
    integer,parameter:: ilist = 45
    integer,parameter:: irdfdat = 46
    integer,parameter:: ixtl = 47
    integer,parameter:: icastep = 48
    integer,parameter:: icmb = 48
    integer,parameter:: imolecule = 49
    integer,parameter:: ires = 50
    integer,parameter:: icmp = 51
    integer,parameter:: icmpout = 52
    integer,parameter:: iylm = 53
    integer,parameter:: iwww = 54
    integer,parameter:: irbonds = 55
    integer,parameter:: irangles = 56
    integer,parameter:: irbout = 57
    integer,parameter:: iraout = 58
    integer,parameter:: istuff = 101           
    integer,parameter:: itmp1 = 201           
    integer,parameter:: itmp2 = 202           
    integer,parameter:: itmp3 = 203           
    integer,parameter:: itmp4 = 204           
    integer,parameter:: itmp5 = 205           
    integer,parameter:: itmp6 = 206           
    integer,parameter:: itmp7 = 207           
    integer,parameter:: itmp8 = 208           
    integer,parameter:: itmp9 = 209           

end module channel_numbers

!===============================================================================

module structure_data

    integer, parameter :: nelements = 111, max_atom_string = 60

    double precision :: a,b,c,alpha,beta,gamma,cell(3,3),scell(3,3),density,approxsize, &
                        rnano(3),rsame,rpdfmax,pdfwidth,zeta,fi,csi,psi,rhollow(3),nanobox(3)
    double precision, allocatable :: xf(:),yf(:),zf(:),xo(:),yo(:),zo(:),spin(:,:),occ(:)
    double precision, allocatable ::xco(:),yco(:),zco(:)
    double precision, allocatable :: xrot(:),yrot(:),zrot(:)!xfrot(:),yfrot(:),zfrot(:)
    double precision, allocatable :: concentration(:)
    double precision :: element_mass(-1:nelements),blength(-1:nelements),rgb(-1:nelements,3), Uiso(-1:nelements)
    double precision :: gridcellsize = 5.0d0
    double precision :: xyzlimits(6)
    double precision :: shiftx,shifty,shiftz
    double precision :: scalenano


    integer, allocatable :: atom_type(:),reference_number(:),reference_cell(:,:),atom_label(:)
    integer, allocatable :: numoftype(:),atom_number(:),gridcell(:,:), &
                            gridatoms(:,:,:,:),neighbours(:,:),norder(:)
    integer, allocatable :: ordertype(:),type_number(:)
    integer :: multiplicity,nsym,natoms,lattice_info(2),ncell(3), &
             ntypes,nspins,mtransform(3,3), ncell_rmc6f(3)
    integer :: n_elements(-1:nelements),n_labels(-1:nelements)

    character(len=max_atom_string), allocatable :: caxyz(:),cspin(:)
    character(len=26), allocatable :: coperators(:)
    character(len=4), allocatable :: atom_name(:)
    character(len=4), allocatable :: aname(:)
    character(len=5), allocatable :: apairs(:)
    character(len=10), allocatable :: catom_label(:)
    character(len=2), allocatable :: element_of_type(:)
    character(len=80), allocatable :: www_neighbours(:)
    character(len=12) :: system
    character(len=1) :: centre
    character(len=2) :: element(-1:nelements),elementuc(-1:nelements),elementlc(-1:nelements)
    character(len=13) :: element_name(-1:nelements)
    character(len=15) :: cshape

    logical :: lcentric
    logical :: input_rmc6f = .false.
    logical, allocatable :: latomuse(:)

end module structure_data

!===============================================================================

module cif_data

    character(len=256), allocatable :: atom_data_lines(:)

end module cif_data

!===============================================================================

module cif_stuff

contains

  subroutine count_cif_items(line,n,natomst)

! Converts the somewhat variable lines containing information about atoms within
! CIF files into one line of data per atom

  use cif_data

  implicit none

  character(len=80) :: line(:),temp
  character(len=20), allocatable :: ctemp(:)
  integer :: n,nlines,i,j,ndata,natomst

  nlines = size(line)

  ndata = 0
  do i = 1,nlines
    temp = trim(adjustl(line(i)))
    do while (len(trim(temp))/=0)
    j = index(temp,' ')
    ndata = ndata + 1
    temp = adjustl(temp(j:))
    end do
  end do

  allocate(ctemp(ndata))

  ndata = 0
  do i = 1,nlines
    temp = trim(adjustl(line(i)))
    do while (len(trim(temp))/=0)
    j = index(temp,' ')
    ndata = ndata + 1
    ctemp(ndata) = trim(adjustl(temp(1:j)))
    temp = adjustl(temp(j:))
    end do
  end do

  natomst = ndata/n
  ndata = 0
  allocate(atom_data_lines(natomst))
  atom_data_lines = ''
  do i = 1,natomst
    do j = 1,n
    ndata = ndata + 1
    atom_data_lines(i) = trim(atom_data_lines(i))//' '//trim(ctemp(ndata))
    end do
    atom_data_lines(i) = adjustl(atom_data_lines(i))
  end do

  if (allocated(ctemp)) deallocate(ctemp)

  return

  end subroutine count_cif_items

end module cif_stuff

!===============================================================================

module arguments

    integer, parameter :: max_arg_length = 20

    integer :: nargs,mbank
    integer, allocatable :: useonly(:)
    character(len=120) :: filename,fileroot,filemag,fieldfile,bondfile,wfolder, &
                          potentialsfile,dwfile,moleculefile,listfile,cmbondfile, &
                          comparefile,alignfile
    character(len=10) :: cinput,oshape
    character(len=max_arg_length), allocatable :: arg(:)
    character(len=300) :: comline,cvacancy_list,creplacement_list,corder_list, cUiso_list
    logical :: lcssr,ldlpoly,lcml,lrmc3,lrmc6f,lrmc6o,lcrystal,lorthog,lrect,lfield
    logical :: lsort,lorder,ldiag,lannotate,lhis,lcrush,lmag,lanal,lsupercell,ltransform, lUiso
    logical :: lexp,lnot,lone,lsize,lnano,lmetadata,lwindows,lunix,lcif,lpdf,lvacancy
    logical :: lreplace, lorder_list, lsupercell_list, lnano_list, lpdfbroaden,lconfigin,lzeta,lUiso_list
    logical :: ldelete, lpotentials, labels, lrelabel, lelabel, luselabels, ladjust_xyz_signs
    logical :: lgulp, lgulpc, latomeye, ldw, lreduce, langles, lrotate, lrotate_list, lhollow, lhollow_list, lorient
    logical :: lcmtx, lreorder, lshift, lmolecules, lxyzlimits, lgasp, llistfile
    logical :: lwritef, lwriteo, lrdfdat, lxtl, lcastep, lcrystalmakerbonds, lcmsort
    logical :: lnanobox_list ! ZYP
    logical :: lshiftx, lshifty, lshiftz, lcompare, lalign, lonecell, lnanoscale
    logical :: lwww, lsilica, lylm, lmolcentre, lorigin, ldlanal, lwww_write, lrmcanal
    double precision :: Uiso_factor


end module arguments

!===============================================================================

module annotations

    character(len=80) :: fmetadata
    character(len=80) :: config_title = ''
    character(len=80) :: metadata_title = ''
    character(len=80) :: metadata_affiliation = ''
    character(len=80) :: metadata_owner = ''
    character(len=80) :: metadata_date = ''
    character(len=80) :: metadata_material = ''
    character(len=80) :: metadata_phase = ''
    character(len=80) :: metadata_formula = ''
    character(len=80) :: metadata_purpose = ''
    character(len=80) :: metadata_keywords = ''
    character(len=80) :: metadata_temperature = ''
    character(len=80) :: metadata_pressure = ''
    character(len=80) :: metadata_note = ''
    character(len=80) :: metadata_comment = ''
    character(len=80) :: metadata_source = ''
    character(len=80) :: metadata_author = ''
    character(len=80) :: metadata_spaceg = ''
    integer :: moves_generated = 0, moves_tried = 0, moves_accepted = 0
    integer :: prior_saves = 0
    logical :: lm_title, lm_owner, lm_material, lm_comment, lm_source

end module annotations

!===============================================================================

module histogram_commons

    integer :: nr,npar,ncoord_0
    double precision ::dr
    integer, allocatable :: histogram(:,:)
    integer, allocatable :: typec_0(:),typen_0(:),rcoord_0(:,:),nneigh(:,:)

end module histogram_commons

!===============================================================================


module utilities

contains

  subroutine convert_coperators(c1,c2,i1,i2,cstring)
! ==================================================

!--------------------------------------------------------------------------
!
! Extracts the information from a string containing the symmetry operations
!
!--------------------------------------------------------------------------

  implicit none

  integer :: i1(3),i2(3),ii1,ii2
  integer :: i,n
  character(len=2) :: c1(3),c2(3),cc1,cc2
  character(len=26) :: buffer,cstring
  ! >>>>>>>>> Yuanpeng -> debugging >>>>>>>>>>>>
  ! The length of the characters defined as >>>>
  ! following is changed from 7 to 15 since >>>>
  ! some symmetry operators containing >>>>>>>>>
  ! fractions are with length longer than 7. >>>
  character(len=15) :: c(3),cc
  ! <<<<<<<<< Yuanpeng finishes here <<<<<<<<<<<

  buffer = adjustl(cstring)
  if (index(buffer,',')>0) then
    n = index(buffer,',')
    c(1) = trim(adjustl(buffer(1:n-1)))
    buffer = adjustl(buffer(n+1:))
    n = index(buffer,',')
    c(2) = trim(adjustl(buffer(1:n-1)))
    buffer = adjustl(buffer(n+1:))
    c(3) = trim(adjustl(buffer))
  else
    n = index(buffer,' ')
    c(1) = trim(adjustl(buffer(1:n-1)))
    buffer = adjustl(buffer(n+1:))
    n = index(buffer,' ')
    c(2) = trim(adjustl(buffer(1:n-1)))
    buffer = adjustl(buffer(n+1:))
    c(3) = trim(adjustl(buffer))
  end if
  do i = 1,3
    cc = c(i)
!...Locate existence of fraction
    n = index(cc,'/')
    if (n>0) then
    read(cc(n-1:n-1),'(i1)') ii1
    read(cc(n+1:n+1),'(i1)') ii2
    if (n==2) then
      cc = cc(4:)
    else
      cc = cc(1:n-2)
    end if
    else
    ii1 = 0 ; ii2 = 1
    end if
    i1(i) = ii1 ; i2(i) = ii2

!...Now look at the x,y,z part
    cc = adjustl(cc)
    n = len(trim(cc))
    if ((cc(n:n)=='+').or.(cc(n:n)=='-')) cc = cc(1:n-1)
    if(len(trim(cc))==1) then
    cc1 = '+'//cc(1:1)
    cc2 = '  '
    end if
    if(len(trim(cc))==2) then
    cc1 = cc(1:2)
    cc2 = '  '
    end if
    if(len(trim(cc))==3) then
    cc1 = '+'//cc(1:1)
    cc2 = cc(2:3)
    end if
    if(len(trim(cc))==4) then
    cc1 = cc(1:2)
    cc2 = cc(3:4)
    end if
    c1(i) = cc1 ; c2(i) = cc2
  end do

  end subroutine convert_coperators



 subroutine operators_to_matrix(c1,c2,ifrac1,ifrac2,symarray)
!============================================================

! Converts text symmetry operators into a matrix 

  implicit none

  integer :: ifrac1(3),ifrac2(3)
  integer :: ix1,ix2,iy1,iy2,iz1,iz22
  character(len=2) :: c1(3),c2(3)
  double precision :: symarray(3,4)

  symarray = 0.0d0

  select case(c1(1))
    case('+x')
      symarray(1,1) = 1.0d0
    case('-x')
      symarray(1,1) = -1.0d0
    case('+y')
      symarray(2,1) = 1.0d0
    case('-y')
      symarray(2,1) = -1.0d0
    case('+z')
      symarray(3,1) = 1.0d0
    case('-z')
      symarray(3,1) = -1.0d0
  end select
  select case(c2(1))
    case('+x')
      symarray(1,1) = 1.0d0
    case('-x')
      symarray(1,1) = -1.0d0
    case('+y')
      symarray(2,1) = 1.0d0
    case('-y')
      symarray(2,1) = -1.0d0
    case('+z')
      symarray(3,1) = 1.0d0
    case('-z')
      symarray(3,1) = -1.0d0
  end select
  select case(c1(2))
    case('+x')
      symarray(1,2) = 1.0d0
    case('-x')
      symarray(1,2) = -1.0d0
    case('+y')
      symarray(2,2) = 1.0d0
    case('-y')
      symarray(2,2) = -1.0d0
    case('+z')
      symarray(3,2) = 1.0d0
    case('-z')
      symarray(3,2) = -1.0d0
  end select
  select case(c2(2))
    case('+x')
      symarray(1,2) = 1.0d0
    case('-x')
      symarray(1,2) = -1.0d0
    case('+y')
      symarray(2,2) = 1.0d0
    case('-y')
      symarray(2,2) = -1.0d0
    case('+z')
      symarray(3,2) = 1.0d0
    case('-z')
      symarray(3,2) = -1.0d0
  end select
  select case(c1(3))
    case('+x')
      symarray(1,3) = 1.0d0
    case('-x')
      symarray(1,3) = -1.0d0
    case('+y')
      symarray(2,3) = 1.0d0
    case('-y')
      symarray(2,3) = -1.0d0
    case('+z')
      symarray(3,3) = 1.0d0
    case('-z')
      symarray(3,3) = -1.0d0
  end select
  select case(c2(3))
    case('+x')
      symarray(1,3) = 1.0d0
    case('-x')
      symarray(1,3) = -1.0d0
    case('+y')
      symarray(2,3) = 1.0d0
    case('-y')
      symarray(2,3) = -1.0d0
    case('+z')
      symarray(3,3) = 1.0d0
    case('-z')
      symarray(3,3) = -1.0d0
  end select
  
  symarray(:,4) = dble(ifrac1)/dble(ifrac2)

 end subroutine operators_to_matrix




   subroutine extract_atom_name_label(stringin,atomname,atomlabel)
!  ===============================================================

!--------------------------------------------------------------------------
!
!  Subroutine to extract atom name and atom label from a single string
!
!--------------------------------------------------------------------------

   integer :: n,nlen,i
   character(len=*) :: stringin,atomname,atomlabel
   
   stringin = adjustl(stringin)

!  Check that the first character is an upper case alphabetic character, and convert to
!  upper case if it is a lower-case character
   if( (ichar(stringin(1:1))>=ichar('a')) .and. (ichar(stringin(1:1))<=ichar('z')) ) then
     stringin(1:1) = char(ichar(stringin(1:1))-ichar('a')+ichar('A'))
   end if
!  Check that the second character is an lower case alphabetic character, and convert to
!  lower case if it is a upper-case character
   if(len(stringin) > 1) then
      if( (ichar(stringin(2:2))>=ichar('A')) .and. (ichar(stringin(2:2))<=ichar('Z')) ) then
         stringin(2:2) = char(ichar(stringin(2:2))-ichar('A')+ichar('a'))
      end if
   end if
  
!  Deduce the element name, as either a 1 or 2 character string
   atomname = ''
   nlen = 1
   if(len(stringin) > 1) then
      if( (ichar(stringin(2:2))>=ichar('a')) .and. (ichar(stringin(2:2))<=ichar('z')) ) nlen = 2
   end if
   atomname = stringin(1:nlen)

   if( (ichar(stringin(1:1))<ichar('A')) .or. (ichar(stringin(1:1))>ichar('Z')) ) then
     write(6,'(a)') 'Element appears to have the wrong characters: '//trim(stringin)
     stop
   end if

!  Check whether there is a label, and return if not
   if ( len_trim(stringin)==nlen ) then
     atomlabel = ''
     return
   end if

!  Deduce the atom label, and stop if the label contains non-numeric characters,
!  except the case where a wildcard label is used
   nlen = nlen + 1
   atomlabel = stringin(nlen:)
   n = len_trim(atomlabel)
   if (atomlabel(1:1)=='*') return
   do i = 1,n
     if (atomlabel(i:i)==' ') cycle
     if ( (ichar(atomlabel(i:i))<ichar('0')) .or. (ichar(atomlabel(i:i))>ichar('9')) ) then
       write(6,'(a)') 'Atom label contains non-numeric characters: '//trim(stringin)
       stop
     end if
   end do
   
   return
   
   end subroutine extract_atom_name_label

double precision function random_normal()

! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.

implicit none

!     Local variables
double precision :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472, &
            r1 = 0.27597, r2 = 0.27846, u, v, x, y, q, half = 0.5d0

!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

do
  call random_number(u)
  call random_number(v)
  v = 1.7156 * (v - half)

!     Evaluate the quadratic form
  x = u - s
  y = ABS(v) - t
  q = x**2 + y*(a*y - b*x)

!     Accept P if inside inner ellipse
  if (q < r1) exit
!     Reject P if outside outer ellipse
  if (q > r2) cycle
!     Reject P if outside acceptance region
  if (v**2 < -4.0*log(u)*u**2) exit
end do

!     Return ratio of P's coordinates as the normal deviate
random_normal = v/u
return

end function random_normal
   
   
! ===============================================================================    
!       Funcitons are borrowed from RMCProfile7
! ===============================================================================

pure function cross_product(a,b)
  implicit none
  ! Use 0-based arrays because it makes the expression simpler
  double precision, intent(in), dimension(0:2) :: a,b
  double precision, dimension(0:2) :: cross_product
  
  integer :: i
  
  do i=0,2
    cross_product(i) = a(mod(i+1,3)) * b(mod(i+2,3)) - a(mod(i+2,3)) * b(mod(i+1,3))
  end do
end function cross_product

pure function matrix33_inverse(mat)
  implicit none
  ! Use 0-based arrays because it makes the expression simpler
  double precision, intent(in) :: mat(0:2,0:2)
  double precision :: matrix33_inverse(0:2,0:2)
  
  double precision :: adj(0:2,0:2)
  double precision :: det
  
  integer :: i
  
  do i=0,2
    adj(i,:) = cross_product(mat(:,mod(i+1,3)), mat(:,mod(i+2,3)))
  end do
  
  det = dot_product(mat(:,1), adj(1,:))
  
  matrix33_inverse = adj / det
end function matrix33_inverse

! ===============================================================================
!   END
! ===============================================================================   

end module utilities

!===============================================================================

module fft_mod

contains

! Fast Fourier/Cosine/Sine Transform
!     dimension   :one
!     data length :power of 2
!     decimation  :frequency
!     radix       :4, 2
!     data        :inplace
!     table       :use
! subroutines
!     cdft: Complex Discrete Fourier Transform
!     rdft: Real Discrete Fourier Transform
!     ddct: Discrete Cosine Transform
!     ddst: Discrete Sine Transform
!     dfct: Cosine Transform of RDFT (Real Symmetric DFT)
!     dfst: Sine Transform of RDFT (Real Anti-symmetric DFT)
!
!
! -------- Complex DFT (Discrete Fourier Transform) --------
!     [definition]
!         <case1>
!             X(k) = sum_j=0^n-1 x(j)*exp(2*pi*i*j*k/n), 0<=k<n
!         <case2>
!             X(k) = sum_j=0^n-1 x(j)*exp(-2*pi*i*j*k/n), 0<=k<n
!         (notes: sum_j=0^n-1 is a summation from j=0 to n-1)
!     [usage]
!         <case1>
!             ip(0) = 0  ! first time only
!             call cdft(2*n, 1, a, ip, w)
!         <case2>
!             ip(0) = 0  ! first time only
!             call cdft(2*n, -1, a, ip, w)
!     [parameters]
!         2*n          :data length (integer)
!                       n >= 1, n = power of 2
!         a(0:2*n-1)   :input/output data (real*8)
!                       input data
!                           a(2*j) = Re(x(j)), 
!                           a(2*j+1) = Im(x(j)), 0<=j<n
!                       output data
!                           a(2*k) = Re(X(k)), 
!                           a(2*k+1) = Im(X(k)), 0<=k<n
!         ip(0:*)      :work area for bit reversal (integer)
!                       length of ip >= 2+sqrt(n)  ; if mod(n,4) = 0
!                                       2+sqrt(n/2); otherwise
!                       ip(0),ip(1) are pointers of the cos/sin table.
!         w(0:n/2-1)   :cos/sin table (real*8)
!                       w(),ip() are initialized if ip(0) = 0.
!     [remark]
!         Inverse of 
!             call cdft(2*n, -1, a, ip, w)
!         is 
!             call cdft(2*n, 1, a, ip, w)
!             do j = 0, 2 * n - 1
!                 a(j) = a(j) / n
!             end do
!         .
!
!
! -------- Real DFT / Inverse of Real DFT --------
!     [definition]
!         <case1> RDFT
!             R(k) = sum_j=0^n-1 a(j)*cos(2*pi*j*k/n), 0<=k<=n/2
!             I(k) = sum_j=0^n-1 a(j)*sin(2*pi*j*k/n), 0<k<n/2
!         <case2> IRDFT (excluding scale)
!             a(k) = R(0)/2 + R(n/2)/2 + 
!                    sum_j=1^n/2-1 R(j)*cos(2*pi*j*k/n) + 
!                    sum_j=1^n/2-1 I(j)*sin(2*pi*j*k/n), 0<=k<n
!     [usage]
!         <case1>
!             ip(0) = 0  ! first time only
!             call rdft(n, 1, a, ip, w)
!         <case2>
!             ip(0) = 0  ! first time only
!             call rdft(n, -1, a, ip, w)
!     [parameters]
!         n            :data length (integer)
!                       n >= 2, n = power of 2
!         a(0:n-1)     :input/output data (real*8)
!                       <case1>
!                           output data
!                               a(2*k) = R(k), 0<=k<n/2
!                               a(2*k+1) = I(k), 0<k<n/2
!                               a(1) = R(n/2)
!                       <case2>
!                           input data
!                               a(2*j) = R(j), 0<=j<n/2
!                               a(2*j+1) = I(j), 0<j<n/2
!                               a(1) = R(n/2)
!         ip(0:*)      :work area for bit reversal (integer)
!                       length of ip >= 2+sqrt(n/2); if mod(n,4) = 2
!                                       2+sqrt(n/4); otherwise
!                       ip(0),ip(1) are pointers of the cos/sin table.
!         w(0:n/2-1)   :cos/sin table (real*8)
!                       w(),ip() are initialized if ip(0) = 0.
!     [remark]
!         Inverse of 
!             call rdft(n, 1, a, ip, w)
!         is 
!             call rdft(n, -1, a, ip, w)
!             do j = 0, n - 1
!                 a(j) = a(j) * 2 / n
!             end do
!         .
!
!
! -------- DCT (Discrete Cosine Transform) / Inverse of DCT --------
!     [definition]
!         <case1> IDCT (excluding scale)
!             C(k) = sum_j=0^n-1 a(j)*cos(pi*j*(k+1/2)/n), 0<=k<n
!         <case2> DCT
!             C(k) = sum_j=0^n-1 a(j)*cos(pi*(j+1/2)*k/n), 0<=k<n
!     [usage]
!         <case1>
!             ip(0) = 0  ! first time only
!             call ddct(n, 1, a, ip, w)
!         <case2>
!             ip(0) = 0  ! first time only
!             call ddct(n, -1, a, ip, w)
!     [parameters]
!         n            :data length (integer)
!                       n >= 2, n = power of 2
!         a(0:n-1)     :input/output data (real*8)
!                       output data
!                           a(k) = C(k), 0<=k<n
!         ip(0:*)      :work area for bit reversal (integer)
!                       length of ip >= 2+sqrt(n/2); if mod(n,4) = 2
!                                       2+sqrt(n/4); otherwise
!                       ip(0),ip(1) are pointers of the cos/sin table.
!         w(0:n*5/4-1) :cos/sin table (real*8)
!                       w(),ip() are initialized if ip(0) = 0.
!     [remark]
!         Inverse of 
!             call ddct(n, -1, a, ip, w)
!         is 
!             a(0) = a(0) / 2
!             call ddct(n, 1, a, ip, w)
!             do j = 0, n - 1
!                 a(j) = a(j) * 2 / n
!             end do
!         .
!
!
! -------- DST (Discrete Sine Transform) / Inverse of DST --------
!     [definition]
!         <case1> IDST (excluding scale)
!             S(k) = sum_j=1^n A(j)*sin(pi*j*(k+1/2)/n), 0<=k<n
!         <case2> DST
!             S(k) = sum_j=0^n-1 a(j)*sin(pi*(j+1/2)*k/n), 0<k<=n
!     [usage]
!         <case1>
!             ip(0) = 0  ! first time only
!             call ddst(n, 1, a, ip, w)
!         <case2>
!             ip(0) = 0  ! first time only
!             call ddst(n, -1, a, ip, w)
!     [parameters]
!         n            :data length (integer)
!                       n >= 2, n = power of 2
!         a(0:n-1)     :input/output data (real*8)
!                       <case1>
!                           input data
!                               a(j) = A(j), 0<j<n
!                               a(0) = A(n)
!                           output data
!                               a(k) = S(k), 0<=k<n
!                       <case2>
!                           output data
!                               a(k) = S(k), 0<k<n
!                               a(0) = S(n)
!         ip(0:*)      :work area for bit reversal (integer)
!                       length of ip >= 2+sqrt(n/2); if mod(n,4) = 2
!                                       2+sqrt(n/4); otherwise
!                       ip(0),ip(1) are pointers of the cos/sin table.
!         w(0:n*5/4-1) :cos/sin table (real*8)
!                       w(),ip() are initialized if ip(0) = 0.
!     [remark]
!         Inverse of 
!             call ddst(n, -1, a, ip, w)
!         is 
!             a(0) = a(0) / 2
!             call ddst(n, 1, a, ip, w)
!             do j = 0, n - 1
!                 a(j) = a(j) * 2 / n
!             end do
!         .
!
!
! -------- Cosine Transform of RDFT (Real Symmetric DFT) --------
!     [definition]
!         C(k) = sum_j=0^n a(j)*cos(pi*j*k/n), 0<=k<=n
!     [usage]
!         ip(0) = 0  ! first time only
!         call dfct(n, a, t, ip, w)
!     [parameters]
!         n            :data length - 1 (integer)
!                       n >= 2, n = power of 2
!         a(0:n)       :input/output data (real*8)
!                       output data
!                           a(k) = C(k), 0<=k<=n
!         t(0:n/2)     :work area (real*8)
!         ip(0:*)      :work area for bit reversal (integer)
!                       length of ip >= 2+sqrt(n/4); if mod(n,4) = 0
!                                       2+sqrt(n/8); otherwise
!                       ip(0),ip(1) are pointers of the cos/sin table.
!         w(0:n*5/8-1) :cos/sin table (real*8)
!                       w(),ip() are initialized if ip(0) = 0.
!     [remark]
!         Inverse of 
!             a(0) = a(0) / 2
!             a(n) = a(n) / 2
!             call dfct(n, a, t, ip, w)
!         is 
!             a(0) = a(0) / 2
!             a(n) = a(n) / 2
!             call dfct(n, a, t, ip, w)
!             do j = 0, n
!                 a(j) = a(j) * 2 / n
!             end do
!         .
!
!
! -------- Sine Transform of RDFT (Real Anti-symmetric DFT) --------
!     [definition]
!         S(k) = sum_j=1^n-1 a(j)*sin(pi*j*k/n), 0<k<n
!     [usage]
!         ip(0) = 0  ! first time only
!         call dfst(n, a, t, ip, w)
!     [parameters]
!         n            :data length + 1 (integer)
!                       n >= 2, n = power of 2
!         a(0:n-1)     :input/output data (real*8)
!                       output data
!                           a(k) = S(k), 0<k<n
!                       (a(0) is used for work area)
!         t(0:n/2-1)   :work area (real*8)
!         ip(0:*)      :work area for bit reversal (integer)
!                       length of ip >= 2+sqrt(n/4); if mod(n,4) = 0
!                                       2+sqrt(n/8); otherwise
!                       ip(0),ip(1) are pointers of the cos/sin table.
!         w(0:n*5/8-1) :cos/sin table (real*8)
!                       w(),ip() are initialized if ip(0) = 0.
!     [remark]
!         Inverse of 
!             call dfst(n, a, t, ip, w)
!         is 
!             call dfst(n, a, t, ip, w)
!             do j = 1, n - 1
!                 a(j) = a(j) * 2 / n
!             end do
!         .
!
!
      subroutine cdft(n, isgn, a, ip, w)
      integer n, isgn, ip(0 : *), j
      real(8) a(0 : n - 1), w(0 : *)
      if (n .gt. 4 * ip(0)) then
          call makewt(n / 4, ip, w)
      end if
      if (n .gt. 4) call bitrv2(n, ip(2), a)
      if (n .gt. 4 .and. isgn .lt. 0) then
          do j = 1, n - 1, 2
              a(j) = -a(j)
          end do
          call cftsub(n, a, w)
          do j = 1, n - 1, 2
              a(j) = -a(j)
          end do
      else
          call cftsub(n, a, w)
      end if
      end subroutine cdft
!
      subroutine rdft(n, isgn, a, ip, w)
      integer n, isgn, ip(0 : *), j, nw, nc
      real(8) a(0 : n - 1), w(0 : *), xi
      nw = ip(0)
      if (n .gt. 4 * nw) then
          nw = n / 4
          call makewt(nw, ip, w)
      end if
      nc = ip(1)
      if (n .gt. 4 * nc) then
          nc = n / 4
          call makect(nc, ip, w(nw))
      end if
      if (isgn .lt. 0) then
          a(1) = 0.5d0 * (a(1) - a(0))
          a(0) = a(0) + a(1)
          do j = 3, n - 1, 2
              a(j) = -a(j)
          end do
          if (n .gt. 4) then
              call rftsub(n, a, nc, w(nw))
              call bitrv2(n, ip(2), a)
          end if
          call cftsub(n, a, w)
          do j = 1, n - 1, 2
              a(j) = -a(j)
          end do
      else
          if (n .gt. 4) call bitrv2(n, ip(2), a)
          call cftsub(n, a, w)
          if (n .gt. 4) call rftsub(n, a, nc, w(nw))
          xi = a(0) - a(1)
          a(0) = a(0) + a(1)
          a(1) = xi
      end if
      end subroutine rdft
!
      subroutine ddct(n, isgn, a, ip, w)
      integer n, isgn, ip(0 : *), j, nw, nc
      real(8) a(0 : n - 1), w(0 : *), xr
      nw = ip(0)
      if (n .gt. 4 * nw) then
          nw = n / 4
          call makewt(nw, ip, w)
      end if
      nc = ip(1)
      if (n .gt. nc) then
          nc = n
          call makect(nc, ip, w(nw))
      end if
      if (isgn .lt. 0) then
          xr = a(n - 1)
          do j = n - 2, 2, -2
              a(j + 1) = a(j - 1) - a(j)
              a(j) = a(j) + a(j - 1)
          end do
          a(1) = xr - a(0)
          a(0) = a(0) + xr
          if (n .gt. 4) then
              call rftsub(n, a, nc, w(nw))
              call bitrv2(n, ip(2), a)
          end if
          call cftsub(n, a, w)
          do j = 1, n - 1, 2
              a(j) = -a(j)
          end do
      end if
      call dctsub(n, a, nc, w(nw))
      if (isgn .ge. 0) then
          if (n .gt. 4) call bitrv2(n, ip(2), a)
          call cftsub(n, a, w)
          if (n .gt. 4) call rftsub(n, a, nc, w(nw))
          xr = a(0) - a(1)
          a(0) = a(0) + a(1)
          do j = 2, n - 2, 2
              a(j - 1) = a(j) - a(j + 1)
              a(j) = a(j) + a(j + 1)
          end do
          a(n - 1) = xr
      end if
      end subroutine ddct
!
      subroutine ddst(n, isgn, a, ip, w)
      integer n, isgn, ip(0 : *), j, nw, nc
      real(8) a(0 : n - 1), w(0 : *), xr
      nw = ip(0)
      if (n .gt. 4 * nw) then
          nw = n / 4
          call makewt(nw, ip, w)
      end if
      nc = ip(1)
      if (n .gt. nc) then
          nc = n
          call makect(nc, ip, w(nw))
      end if
      if (isgn .lt. 0) then
          xr = a(n - 1)
          do j = n - 2, 2, -2
              a(j + 1) = a(j - 1) + a(j)
              a(j) = a(j) - a(j - 1)
          end do
          a(1) = -xr - a(0)
          a(0) = a(0) - xr
          if (n .gt. 4) then
              call rftsub(n, a, nc, w(nw))
              call bitrv2(n, ip(2), a)
          end if
          call cftsub(n, a, w)
          do j = 1, n - 1, 2
              a(j) = -a(j)
          end do
      end if
      call dstsub(n, a, nc, w(nw))
      if (isgn .ge. 0) then
          if (n .gt. 4) call bitrv2(n, ip(2), a)
          call cftsub(n, a, w)
          if (n .gt. 4) call rftsub(n, a, nc, w(nw))
          xr = a(0) - a(1)
          a(0) = a(0) + a(1)
          do j = 2, n - 2, 2
              a(j - 1) = -a(j) - a(j + 1)
              a(j) = a(j) - a(j + 1)
          end do
          a(n - 1) = -xr
      end if
      end subroutine ddst
!
      subroutine dfct(n, a, t, ip, w)
      integer n, ip(0 : *), j, k, l, m, mh, nw, nc
      real(8) a(0 : n), t(0 : n / 2), w(0 : *), xr, xi
      nw = ip(0)
      if (n .gt. 8 * nw) then
          nw = n / 8
          call makewt(nw, ip, w)
      end if
      nc = ip(1)
      if (n .gt. 2 * nc) then
          nc = n / 2
          call makect(nc, ip, w(nw))
      end if
      m = n / 2
      xr = a(0) + a(n)
      a(0) = a(0) - a(n)
      t(0) = xr - a(m)
      t(m) = xr + a(m)
      if (n .gt. 2) then
          mh = m / 2
          do j = 1, mh - 1
              k = m - j
              xr = a(j) + a(n - j)
              a(j) = a(j) - a(n - j)
              xi = a(k) + a(n - k)
              a(k) = a(k) - a(n - k)
              t(j) = xr - xi
              t(k) = xr + xi
          end do
          t(mh) = a(mh) + a(n - mh)
          a(mh) = a(mh) - a(n - mh)
          call dctsub(m, a, nc, w(nw))
          if (m .gt. 4) call bitrv2(m, ip(2), a)
          call cftsub(m, a, w)
          if (m .gt. 4) call rftsub(m, a, nc, w(nw))
          xr = a(0) + a(1)
          a(n - 1) = a(0) - a(1)
          do j = m - 2, 2, -2
              a(2 * j + 1) = a(j) + a(j + 1)
              a(2 * j - 1) = a(j) - a(j + 1)
          end do
          a(1) = xr
          l = 2
          m = mh
          do while (m .ge. 2)
              call dctsub(m, t, nc, w(nw))
              if (m .gt. 4) call bitrv2(m, ip(2), t)
              call cftsub(m, t, w)
              if (m .gt. 4) call rftsub(m, t, nc, w(nw))
              a(n - l) = t(0) - t(1)
              a(l) = t(0) + t(1)
              k = 0
              do j = 2, m - 2, 2
                  k = k + 4 * l
                  a(k - l) = t(j) - t(j + 1)
                  a(k + l) = t(j) + t(j + 1)
              end do
              l = 2 * l
              mh = m / 2
              do j = 0, mh - 1
                  k = m - j
                  t(j) = t(m + k) - t(m + j)
                  t(k) = t(m + k) + t(m + j)
              end do
              t(mh) = t(m + mh)
              m = mh
          end do
          a(l) = t(0)
          a(n) = t(2) - t(1)
          a(0) = t(2) + t(1)
      else
          a(1) = a(0)
          a(2) = t(0)
          a(0) = t(1)
      end if
      end subroutine dfct
!
      subroutine dfst(n, a, t, ip, w)
      integer n, ip(0 : *), j, k, l, m, mh, nw, nc
      real(8) a(0 : n - 1), t(0 : n / 2 - 1), w(0 : *), xr, xi
      nw = ip(0)
      if (n .gt. 8 * nw) then
          nw = n / 8
          call makewt(nw, ip, w)
      end if
      nc = ip(1)
      if (n .gt. 2 * nc) then
          nc = n / 2
          call makect(nc, ip, w(nw))
      end if
      if (n .gt. 2) then
          m = n / 2
          mh = m / 2
          do j = 1, mh - 1
              k = m - j
              xr = a(j) - a(n - j)
              a(j) = a(j) + a(n - j)
              xi = a(k) - a(n - k)
              a(k) = a(k) + a(n - k)
              t(j) = xr + xi
              t(k) = xr - xi
          end do
          t(0) = a(mh) - a(n - mh)
          a(mh) = a(mh) + a(n - mh)
          a(0) = a(m)
          call dstsub(m, a, nc, w(nw))
          if (m .gt. 4) call bitrv2(m, ip(2), a)
          call cftsub(m, a, w)
          if (m .gt. 4) call rftsub(m, a, nc, w(nw))
          xr = a(0) + a(1)
          a(n - 1) = a(1) - a(0)
          do j = m - 2, 2, -2
              a(2 * j + 1) = a(j) - a(j + 1)
              a(2 * j - 1) = -a(j) - a(j + 1)
          end do
          a(1) = xr
          l = 2
          m = mh
          do while (m .ge. 2)
              call dstsub(m, t, nc, w(nw))
              if (m .gt. 4) call bitrv2(m, ip(2), t)
              call cftsub(m, t, w)
              if (m .gt. 4) call rftsub(m, t, nc, w(nw))
              a(n - l) = t(1) - t(0)
              a(l) = t(0) + t(1)
              k = 0
              do j = 2, m - 2, 2
                  k = k + 4 * l
                  a(k - l) = -t(j) - t(j + 1)
                  a(k + l) = t(j) - t(j + 1)
              end do
              l = 2 * l
              mh = m / 2
              do j = 1, mh - 1
                  k = m - j
                  t(j) = t(m + k) + t(m + j)
                  t(k) = t(m + k) - t(m + j)
              end do
              t(0) = t(m + mh)
              m = mh
          end do
          a(l) = t(0)
      end if
      a(0) = 0
      end subroutine dfst
!
! -------- initializing routines --------
!
      subroutine makewt(nw, ip, w)
      integer nw, ip(0 : *), nwh, j
      real(8) w(0 : nw - 1), delta, x, y
      ip(0) = nw
      ip(1) = 1
      if (nw .gt. 2) then
          nwh = nw / 2
          delta = atan(1.0d0) / nwh
          w(0) = 1
          w(1) = 0
          w(nwh) = cos(delta * nwh)
          w(nwh + 1) = w(nwh)
          do j = 2, nwh - 2, 2
              x = cos(delta * j)
              y = sin(delta * j)
              w(j) = x
              w(j + 1) = y
              w(nw - j) = y
              w(nw - j + 1) = x
          end do
          call bitrv2(nw, ip(2), w)
      end if
      end subroutine makewt
!
      subroutine makect(nc, ip, c)
      integer nc, ip(0 : *), nch, j
      real(8) c(0 : nc - 1), delta
      ip(1) = nc
      if (nc .gt. 1) then
          nch = nc / 2
          delta = atan(1.0d0) / nch
          c(0) = 0.5d0
          c(nch) = 0.5d0 * cos(delta * nch)
          do j = 1, nch - 1
              c(j) = 0.5d0 * cos(delta * j)
              c(nc - j) = 0.5d0 * sin(delta * j)
          end do
      end if
      end subroutine makect
!
! -------- child routines --------
!
      subroutine bitrv2(n, ip, a)
      integer n, ip(0 : *), j, j1, k, k1, l, m, m2
      real(8) a(0 : n - 1), xr, xi
      ip(0) = 0
      l = n
      m = 1
      do while (4 * m .lt. l)
          l = l / 2
          do j = 0, m - 1
              ip(m + j) = ip(j) + l
          end do
          m = m * 2
      end do
      if (4 * m .gt. l) then
          do k = 1, m - 1
              do j = 0, k - 1
                  j1 = 2 * j + ip(k)
                  k1 = 2 * k + ip(j)
                  xr = a(j1)
                  xi = a(j1 + 1)
                  a(j1) = a(k1)
                  a(j1 + 1) = a(k1 + 1)
                  a(k1) = xr
                  a(k1 + 1) = xi
              end do
          end do
      else
          m2 = 2 * m
          do k = 1, m - 1
              do j = 0, k - 1
                  j1 = 2 * j + ip(k)
                  k1 = 2 * k + ip(j)
                  xr = a(j1)
                  xi = a(j1 + 1)
                  a(j1) = a(k1)
                  a(j1 + 1) = a(k1 + 1)
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + m2
                  k1 = k1 + m2
                  xr = a(j1)
                  xi = a(j1 + 1)
                  a(j1) = a(k1)
                  a(j1 + 1) = a(k1 + 1)
                  a(k1) = xr
                  a(k1 + 1) = xi
              end do
          end do
      end if
      end subroutine bitrv2
!
      subroutine cftsub(n, a, w)
      integer n, j, j1, j2, j3, k, k1, ks, l, m
      real(8) a(0 : n - 1), w(0 : *)
      real(8) wk1r, wk1i, wk2r, wk2i, wk3r, wk3i
      real(8) x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
      l = 2
      do while (2 * l .lt. n)
          m = 4 * l
          do j = 0, l - 2, 2
              j1 = j + l
              j2 = j1 + l
              j3 = j2 + l
              x0r = a(j) + a(j1)
              x0i = a(j + 1) + a(j1 + 1)
              x1r = a(j) - a(j1)
              x1i = a(j + 1) - a(j1 + 1)
              x2r = a(j2) + a(j3)
              x2i = a(j2 + 1) + a(j3 + 1)
              x3r = a(j2) - a(j3)
              x3i = a(j2 + 1) - a(j3 + 1)
              a(j) = x0r + x2r
              a(j + 1) = x0i + x2i
              a(j2) = x0r - x2r
              a(j2 + 1) = x0i - x2i
              a(j1) = x1r - x3i
              a(j1 + 1) = x1i + x3r
              a(j3) = x1r + x3i
              a(j3 + 1) = x1i - x3r
          end do
          if (m .lt. n) then
              wk1r = w(2)
              do j = m, l + m - 2, 2
                  j1 = j + l
                  j2 = j1 + l
                  j3 = j2 + l
                  x0r = a(j) + a(j1)
                  x0i = a(j + 1) + a(j1 + 1)
                  x1r = a(j) - a(j1)
                  x1i = a(j + 1) - a(j1 + 1)
                  x2r = a(j2) + a(j3)
                  x2i = a(j2 + 1) + a(j3 + 1)
                  x3r = a(j2) - a(j3)
                  x3i = a(j2 + 1) - a(j3 + 1)
                  a(j) = x0r + x2r
                  a(j + 1) = x0i + x2i
                  a(j2) = x2i - x0i
                  a(j2 + 1) = x0r - x2r
                  x0r = x1r - x3i
                  x0i = x1i + x3r
                  a(j1) = wk1r * (x0r - x0i)
                  a(j1 + 1) = wk1r * (x0r + x0i)
                  x0r = x3i + x1r
                  x0i = x3r - x1i
                  a(j3) = wk1r * (x0i - x0r)
                  a(j3 + 1) = wk1r * (x0i + x0r)
              end do
              k1 = 1
              ks = -1
              do k = 2 * m, n - m, m
                  k1 = k1 + 1
                  ks = -ks
                  wk1r = w(2 * k1)
                  wk1i = w(2 * k1 + 1)
                  wk2r = ks * w(k1)
                  wk2i = w(k1 + ks)
                  wk3r = wk1r - 2 * wk2i * wk1i
                  wk3i = 2 * wk2i * wk1r - wk1i
                  do j = k, l + k - 2, 2
                      j1 = j + l
                      j2 = j1 + l
                      j3 = j2 + l
                      x0r = a(j) + a(j1)
                      x0i = a(j + 1) + a(j1 + 1)
                      x1r = a(j) - a(j1)
                      x1i = a(j + 1) - a(j1 + 1)
                      x2r = a(j2) + a(j3)
                      x2i = a(j2 + 1) + a(j3 + 1)
                      x3r = a(j2) - a(j3)
                      x3i = a(j2 + 1) - a(j3 + 1)
                      a(j) = x0r + x2r
                      a(j + 1) = x0i + x2i
                      x0r = x0r - x2r
                      x0i = x0i - x2i
                      a(j2) = wk2r * x0r - wk2i * x0i
                      a(j2 + 1) = wk2r * x0i + wk2i * x0r
                      x0r = x1r - x3i
                      x0i = x1i + x3r
                      a(j1) = wk1r * x0r - wk1i * x0i
                      a(j1 + 1) = wk1r * x0i + wk1i * x0r
                      x0r = x1r + x3i
                      x0i = x1i - x3r
                      a(j3) = wk3r * x0r - wk3i * x0i
                      a(j3 + 1) = wk3r * x0i + wk3i * x0r
                  end do
              end do
          end if
          l = m
      end do
      if (l .lt. n) then
          do j = 0, l - 2, 2
              j1 = j + l
              x0r = a(j) - a(j1)
              x0i = a(j + 1) - a(j1 + 1)
              a(j) = a(j) + a(j1)
              a(j + 1) = a(j + 1) + a(j1 + 1)
              a(j1) = x0r
              a(j1 + 1) = x0i
          end do
      end if
      end subroutine cftsub
!
      subroutine rftsub(n, a, nc, c)
      integer n, nc, j, k, kk, ks
      real(8) a(0 : n - 1), c(0 : nc - 1), wkr, wki, xr, xi, yr, yi
      ks = 4 * nc / n
      kk = 0
      do k = n / 2 - 2, 2, -2
          j = n - k
          kk = kk + ks
          wkr = 0.5d0 - c(kk)
          wki = c(nc - kk)
          xr = a(k) - a(j)
          xi = a(k + 1) + a(j + 1)
          yr = wkr * xr - wki * xi
          yi = wkr * xi + wki * xr
          a(k) = a(k) - yr
          a(k + 1) = a(k + 1) - yi
          a(j) = a(j) + yr
          a(j + 1) = a(j + 1) - yi
      end do
      end subroutine rftsub
!
      subroutine dctsub(n, a, nc, c)
      integer n, nc, j, k, kk, ks, m
      real(8) a(0 : n - 1), c(0 : nc - 1), wkr, wki, xr
      ks = nc / n
      kk = ks
      m = n / 2
      do k = 1, m - 1
          j = n - k
          wkr = c(kk) - c(nc - kk)
          wki = c(kk) + c(nc - kk)
          kk = kk + ks
          xr = wki * a(k) - wkr * a(j)
          a(k) = wkr * a(k) + wki * a(j)
          a(j) = xr
      end do
      a(m) = 2 * c(kk) * a(m)
      end subroutine dctsub
!
      subroutine dstsub(n, a, nc, c)
      integer n, nc, j, k, kk, ks, m
      real(8) a(0 : n - 1), c(0 : nc - 1), wkr, wki, xr
      ks = nc / n
      kk = ks
      m = n / 2
      do k = 1, m - 1
          j = n - k
          wkr = c(kk) - c(nc - kk)
          wki = c(kk) + c(nc - kk)
          kk = kk + ks
          xr = wki * a(j) - wkr * a(k)
          a(j) = wkr * a(j) + wki * a(k)
          a(k) = xr
      end do
      a(m) = 2 * c(kk) * a(m)
      end subroutine dstsub
!
    end module fft_mod
 
    
    
      BLOCK DATA ZBQLBD01
!*
!*       Initializes seed array etc. for random number generator.
!*       The values below have themselves been generated using the
!*       NAG generator.
!*
      COMMON /ZBQL0001/ ZBQLIX,B,C
      DOUBLE PRECISION ZBQLIX(43),B,C
      INTEGER I
      DATA (ZBQLIX(I),I=1,43) /8.001441D7,5.5321801D8,       &
      1.69570999D8,2.88589940D8,2.91581871D8,1.03842493D8,   &
      7.9952507D7,3.81202335D8,3.11575334D8,4.02878631D8,    &
      2.49757109D8,1.15192595D8,2.10629619D8,3.99952890D8,   &
      4.12280521D8,1.33873288D8,7.1345525D7,2.23467704D8,    &
      2.82934796D8,9.9756750D7,1.68564303D8,2.86817366D8,    &
      1.14310713D8,3.47045253D8,9.3762426D7 ,1.09670477D8,   &
      3.20029657D8,3.26369301D8,9.441177D6,3.53244738D8,     &
      2.44771580D8,1.59804337D8,2.07319904D8,3.37342907D8,   &
      3.75423178D8,7.0893571D7 ,4.26059785D8,3.95854390D8,   &
      2.0081010D7,5.9250059D7,1.62176640D8,3.20429173D8,     &
      2.63576576D8/
      DATA B / 4.294967291D9 /
      DATA C / 0.0D0 /
      END
!******************************************************************
!******************************************************************
!******************************************************************    
    
    
module random
 
! The functions listed below are taken from external source (indicated below) and adopted to be used here. 
! In fact the function ZBQLNOR is the rangom generator from Gaussian distribution.
! It s going to be used for -Uiso option to randomise atom positions from the average accordieng to the given Uiso values.
  
!*******************************************************************
!********	FILE: randgen.f			                     ***********
!********	AUTHORS: Richard Chandler		             ***********
!********	(richard@stats.ucl.ac.uk)	                 ***********
!********	Paul Northrop 			                     ***********
!********	(northrop@stats.ox.ac.uk)	                 ***********
!********	LAST MODIFIED: 26/8/03			             ***********
!********	See file randgen.txt for details             ***********
!*******************************************************************  
!
!    http://www.ucl.ac.uk/~ucakarc/work/randgen.html

    contains 
    
    
    
    !******************************************************************
      FUNCTION ZBQLU01(DUMMY)
!*
!*       Returns a uniform random number between 0 & 1, using
!*       a Marsaglia-Zaman type subtract-with-borrow generator.
!*       Uses double precision, rather than integer, arithmetic 
!*       throughout because MZ's integer constants overflow
!*       32-bit integer storage (which goes from -2^31 to 2^31).
!*       Ideally, we would explicitly truncate all integer 
!*       quantities at each stage to ensure that the double
!*       precision representations do not accumulate approximation
!*       error; however, on some machines the use of DNINT to
!*       accomplish this is *seriously* slow (run-time increased
!*       by a factor of about 3). This double precision version 
!*       has been tested against an integer implementation that
!*       uses long integers (non-standard and, again, slow) -
!*       the output was identical up to the 16th decimal place
!*       after 10^10 calls, so we're probably OK ...
!*
      DOUBLE PRECISION ZBQLU01,DUMMY,B,C,ZBQLIX(43),X,B2,BINV
      INTEGER CURPOS,ID22,ID43

      COMMON /ZBQL0001/ ZBQLIX,B,C
      SAVE /ZBQL0001/
      SAVE CURPOS,ID22,ID43
      DATA CURPOS,ID22,ID43 /1,22,43/

      B2 = B
      BINV = 1.0D0/B
 5    X = ZBQLIX(ID22) - ZBQLIX(ID43) - C
      IF (X.LT.0.0D0) THEN
       X = X + B
       C = 1.0D0
      ELSE
       C = 0.0D0
      ENDIF
      ZBQLIX(ID43) = X
!*
!*     Update array pointers. Do explicit check for bounds of each to
!*     avoid expense of modular arithmetic. If one of them is 0 the others
!*     won't be
!*
      CURPOS = CURPOS - 1
      ID22 = ID22 - 1
      ID43 = ID43 - 1
      IF (CURPOS.EQ.0) THEN
       CURPOS=43
      ELSEIF (ID22.EQ.0) THEN
       ID22 = 43
      ELSEIF (ID43.EQ.0) THEN
       ID43 = 43
      ENDIF
!*
!*     The integer arithmetic there can yield X=0, which can cause 
!*     problems in subsequent routines (e.g. ZBQLEXP). The problem
!*     is simply that X is discrete whereas U is supposed to 
!*     be continuous - hence if X is 0, go back and generate another
!*     X and return X/B^2 (etc.), which will be uniform on (0,1/B). 
!*
      IF (X.LT.BINV) THEN
       B2 = B2*B
       GOTO 5
      ENDIF

      ZBQLU01 = X/B2

      END FUNCTION ZBQLU01 
      

!  ******************************************************************
   FUNCTION ZBQLNOR(MU,SIGMA)
!*
!*       Returns a random number Normally distributed with mean
!*       MU and standard deviation |SIGMA|, using the Box-Muller
!*       algorithm
!*
      
      DOUBLE PRECISION THETA,R,ZBQLNOR,PI,MU,SIGMA
      DOUBLE PRECISION SPARE
      INTEGER STATUS
      SAVE STATUS,SPARE,PI
      DATA STATUS /-1/

      IF (STATUS.EQ.-1) PI = 4.0D0*DATAN(1.0D0)

      IF (STATUS.LE.0) THEN
       THETA = 2.0D0*PI*ZBQLU01(0.0D0)
       R = DSQRT( -2.0D0*DLOG(ZBQLU01(0.0D0)) )
       ZBQLNOR = (R*DCOS(THETA))
       SPARE = (R*DSIN(THETA))
       STATUS = 1
      ELSE
       ZBQLNOR = SPARE
       STATUS = 0
      ENDIF
      
      ZBQLNOR = MU + (SIGMA*ZBQLNOR)

      END FUNCTION ZBQLNOR
      
end module random    
    
    
    

!===============================================================================

  program data2config
! ===================

    use structure_data
    use arguments
    use version_data
    use channel_numbers

    implicit none

    integer :: i,n,ios
    character(len=132) :: buffer

    wfolder = ''

    call get_arguments

    ios = 0 ; open(ic,file=trim(filename),form='formatted',status='old',iostat=ios)
    if (ios/=0) then
      write(6,'(a)') 'Error opening data file'
      stop
    end if

    n = index(filename,'.')
    if ((trim(filename(n+1:))=='cfg').or.(trim(filename(n+1:))=='CFG') ) then
!     Check which type of .cfg file we have; RMCprofile or AtomEye. The only way is to read
!     the first line
      read(ic,'(a)') buffer
      rewind(ic)
      if (buffer(1:21)=='Number of particles =') then
        cinput = 'atomeye'
        latomeye = .false.
        write(6,'(a)') trim(filename)//' is identified as an atomeye file'
      else
        cinput = 'cfg'
        lrmc3 = .false.
        write(6,'(a)') trim(filename)//' is identified as a classic RMC file'
      end if
    end if

    natoms = 0
    ntypes = 0
    nspins = 0
    n_elements = 0
    n_labels = 0

    if (ldiag) then
      n = index(fileroot,'.')
      if (n>0) then
        open(main_output,file=trim(wfolder)//fileroot(1:n)//'out',form='formatted',status='unknown')
      else
        open(main_output,file=trim(wfolder)//trim(fileroot)//'.out',form='formatted',status='unknown')
      end if
      inquire(main_output,name=buffer)
      write(6,'(a)') 'Diagfile = '//trim(buffer)
      write(main_output,'(a/)') 'data2config, version '//trim(version_number(nversions)) &
           //' ('//trim(version_date(nversions))//')'
    end if

    call assign_element_names
    if (lUiso) call assign_element_Uiso
    if (lcssr.or.ldlpoly.or.lrmc6f.or.lrmc6o) call annotate

    if (cinput(1:3)=='tbl') then
      call get_symmetry_tbl
      if (lmag) then
        call get_exp_spins
      else
        call get_structure_tbl
      end if
      call label_atoms
      call generate_structure
      if (lexp.and.(.not.lnot)) call get_expfile_parameters
    end if

    if (cinput(1:3)=='cif') then
      call get_structure_cif
      call label_atoms
      call generate_structure
    end if

    if (cinput(1:3)=='www') call obtain_structure_www

    if (cinput(1:3)=='cfg') then
      if (lmag) call get_cfg_spins
      call obtain_structure_cfg
    end if

    if (cinput(1:5)=='rmc6f') then
!    if (lmag) call get_cfg_spins
      input_rmc6f = .true.
      call obtain_structure_rmc6f
    end if

    if (cinput(1:3)=='his') call obtain_structure_his

    if (cinput(1:4)=='cssr') call obtain_structure_cssr

    if (cinput(1:4)=='cell') call obtain_structure_cell
    
    if (cinput(1:3)=='sff') then
      call get_structure_sff
      call label_atoms
      call generate_structure
    end if

    if (cinput(1:7)=='crystal') call obtain_structure_crystal

    if (cinput(1:3)=='res') then
      call get_structure_resfile
      call label_atoms
      call generate_structure
    end if
    
    if ((cinput(1:6)=='CONFIG').or.(cinput(1:6)=='REVCON')) then
      call obtain_structure_dlpoly
      inquire(file='RDFDAT',exist=lrdfdat)
    end if

    if (trim(cinput)=='atomeye') call obtain_structure_atomeye

    close(ic)

    if (ldiag) then
    write(main_output,*)
    write(main_output,'(/a)') 'Output from data2config'
    write(main_output,'(a/)') '======================='
    if (allocated(coperators)) then
      do i = 1,nsym ; write(main_output,'(a)') coperators(i) ; end do
    end if
    if (cinput(1:3)=='tbl') write(main_output,'(a)') centre
    if (cinput(1:3)=='tbl') write(main_output,*) lcentric
    if (cinput(1:3)=='tbl') write(main_output,*) system
    write(main_output,*) a,b,c
    write(main_output,*) alpha,beta,gamma
!      write(main_output,*) natoms,shape(caxyz)
!      do i = 1,natoms ; write(main_output,'(a)') caxyz(i) ; end do
    end if

    if (lrect) call rectangulate  

    size: if (lsize) then
      ncell(1) = nint(approxsize/a)
      ncell(2) = nint(approxsize/b)
      ncell(3) = nint(approxsize/c)
!       ncell(1) = nint(approxsize/a)*2 ! ZYP -> The original approximate size for the supercell 
!       ncell(2) = nint(approxsize/b)*2 ! ZYP -> sometimes may not be big enough, resulting in 
!       ncell(3) = nint(approxsize/c)*2 ! ZYP -> weird truncation shape of the generated nanoparticle.
    else if ((.not.lsupercell_list).and.(.not.lone).and.(.not.lnano).and.(.not.lreorder)) then
      write(6,'(a)', advance="no") 'Please give supercell multipliers: '
      read(5,*) ncell
    end if size

    nanosize: if (lnano.and.(.not.lnano_list)) then
      write(6,'(a)',advance="no") 'Please give nanoparticle dimensions: '
      read(5,*) rnano 
    end if nanosize
      
    if (lhollow.and..not.lhollow_list) then
      write(6,'(a)') 'Nanoparticle hole dimensions have to be 3 positive numbers!'
      write(6,'(a)',advance="no") 'Please provide them here:'
      read(5,*) rhollow
    end if
      
    if (lrotate.and..not.lrotate_list) then
      write(6,'(a)',advance="no") 'Please provide three Euler angles of rotation: '
      read(5,*) fi,csi,psi
    end if
    
    if (ldelete) lsort = .true.
    if (lshiftx) lshift = .true.
    if (lshifty) lshift = .true.
    if (lshiftz) lshift = .true.
    if (lorigin) call reset_origin

    if (lwww.and.lsilica) call make_into_silica
    if (.not.lone) call expand_structure !expand unit cell
    if (lrotate) call rotate_structure
    if (lnano) call nanoparticle !and then chop it
    if (ltransform) call transform_structure
    if (lreplace) call replace_atoms
    if (lvacancy) call create_vacancies
    if (lsort) call sort_structure
    if (lorder.or.lreorder) call order_structure
    if (lreduce) call reduce
    if (lshift) call shift_origin
    if (lalign) call align_files
    call allocate_grid
    if (lmolecules) call locate_molecules
    if (lorient) call orient_molecules
    if (lnanoscale) call scale_nanoparticle
    call allocate_grid

    if (lUiso) call Uiso_structure
    call statistics 


    if (ldiag) then
      write(main_output,'(/a)') 'Output before call to write routines'
      write(main_output,'(a/)') '===================================='
      write(main_output,'(/a/)') 'Contents of cxyz'
      write(main_output,'(/a/)') 'Element, label, f and o coordinates'
      write(main_output,*) a,b,c
      write(main_output,*) alpha,beta,gamma
      do i = 1,natoms
        write(main_output,'(a,6f12.4)') element(atom_type(i)),xf(i),yf(i),zf(i),xo(i),yo(i),zo(i)
      end do
    end if

    if (cinput(1:7)=='HISTORY') then
      call obtain_structures_history
      stop
    end if

    allocate(latomuse(natoms))
    latomuse = .true.
    if (lxyzlimits) call assign_atoms_to_use

    if (lcssr) call write_cssr
    if (lcmtx) call write_cmtx
    if (lcif) call write_cif
    if (lxtl) call write_xtl
    if (ldlpoly) call write_dlpoly
    if (ldw) call write_distance_windows_file
    if (lgulp) call write_gulp
    if (latomeye) call write_atomeye
    if (lgasp) call write_gasp
    if (lrmc3) call write_rmc3
    if (lrmc6f) call write_rmc6f
    if (lrmc6o) call write_rmc6o
    if (lwww_write) call write_www
    if (lhis) then
      call write_his6f
      call write_his6o
    end if
    if (lcompare) call compare_files
    if (lcrystal) call write_crystal
    if (lcrush) call write_crush
    if (lanal) call analysis
    if (ldlanal) call dlpoly_analysis
    if (lrmcanal) call rmc_analysis
    if (lpdf) call calculate_pdf
    if (lrdfdat) call dlpoly_pdf
    if (lcastep) call write_castep_cell
    if (lcrystalmakerbonds) call crystalmakerbonds
    if (lylm) call ylm

  end program data2config

!===============================================================================

  subroutine get_arguments
! ========================

!--------------------------------------------------------------------------
!
! Reads the arguments when the program is executed as a shell command
!
!--------------------------------------------------------------------------

    use arguments
    use annotations
    use structure_data
    use version_data
    use channel_numbers

    implicit none

    integer :: n,i,ierror,n1,n2,n3,iusehist,nusehist
    character(len=80) :: filetest,clength,ctemp,ctempcount,buffer
    character(len=120) :: temp
    logical :: lokay,l1,l2,lexist
    double precision :: pi

    nargs = command_argument_count()
    allocate(arg(nargs))
    call get_command(comline)

    filename = 'no filename given'
    lcssr = .false. ; ldlpoly = .false. ; lcml = .false. ; lrmc3 = .false.
    lrmc6f = .false. ; lrmc6o = .false. ; lcrystal = .false. ; lrect = .false.
    lsort = .false. ; ldiag = .false. ; lannotate = .true. ; lcrush = .false.
    lmag = .false. ; lorder = .false. ; lexp = .false. ; lnot = .false.
    lone = .false. ; lsize = .false. ; lnano = .false. ; lnano_list = .false. ; lmetadata = .false.
    lsupercell = .false. ; lwindows = .false. ; lunix = .true. ; lcif = .false.
    lpdf = .false. ; lvacancy = .false. ; lorder_list = .false. ; lsupercell_list = .false.
    lpdfbroaden = .false. ; lconfigin = .false. ; lzeta = .false. ; ldelete = .false. ;
    labels = .false. ; lrelabel = .false. ; luselabels = .false. ; ladjust_xyz_signs = .true.
    lgulp = .false. ; lgulpc = .false. ; latomeye = .false. ; ldw = .false. ; lreduce = .false.
    langles = .false. ; lrotate = .false. ; lrotate_list = .false. ; lhollow_list = .false.
    lcmtx = .false. ; lreorder = .false. ; lshift = .false. ; lmolecules = .false. ; lorient = .false.
    lxyzlimits = .false. ; lgasp = .false. ; llistfile = .false.
    lwritef = .true. ; lwriteo = .false. ; lrdfdat = .false. ; lxtl = .false.
    lcastep = .false. ; lcrystalmakerbonds = .false. ; lcmsort = .false.
    lshiftx = .false. ; lshifty = .false. ; lshiftz = .false. ; lcompare = .false.
    lalign = .false. ; lonecell = .false. ; lnanoscale = .false.
    lwww = .false. ; lsilica = .false. ; lorigin = .false. ; lwww_write = .false.
    lUiso_list = .false. ; lUiso = .false. ; lylm = .false. ; ldlanal = .false.
    lrmcanal = .false.

!   Parameter to determine whether two atoms are in the same place; this
!   can be changed using the -rsame parameter set here
    rsame = 0.01d0
    pdfwidth = 0.4d0
    pi = acos(-1.0d0)
!   Default nanoparticle shape 
    cshape = 'sphere' 
    
    mbank = 0

    if (nargs==0) then
    write(6,'(a/)') 'data2config, version '//trim(version_number(nversions)) &
           //' ('//trim(version_date(nversions))//')'
    write(6,*)
    write(6,'(a)') 'Usage: data2config <optional arguments> filename'
    write(6,'(a)') '------------------------'
    write(6,'(a)') 'Optional arguments are:'
    write(6,'(a)') '-help        request for detailed help instructions'
!    write(6,'(a)') '-cml         write CML file'
!    write(6,'(a)') '-crush       write crush input file'
    write(6,'(a)') '-analysis    generate simple bond analysis, requires file input'
    write(6,'(a)') '-angles      generate simple bond angle analysis when -angles selected'
    write(6,'(a)') '-atomeye     write configuration file in Atomeye format'
    write(6,'(a)') '-bank        provide the GSAS histogram number as argument'
    write(6,'(a)') '-castep      write CASTEP cell file'
    write(6,'(a)') '-cif         write cif output file'
    write(6,'(a)') '-cmbonds     bondlength analysis using CrystalMaker distances file, specifed within [...] brackets'
    write(6,'(a)') '-cmsort      sort atoms according to bonds identified by CrystalMaker'
    write(6,'(a)') '-cmtx        write CrystalMaker text input file'
    write(6,'(a)') '-compare     compare the input file with a given file as argument, differences in position given in output'
    write(6,'(a)') '-crystal     write crystal input file'
    write(6,'(a)') '-cssr        write cssr output file'
    write(6,'(a)') '-delete      delete vacant sites if provided in the structure file'
    write(6,'(a)') '-diag        write diagnostics file'
    write(6,'(a)') '-dlanal      Analyse DLPOLY configuration file using the FIELD file'
    write(6,'(a)') '-dlpoly      write DLPOLY CONFIG file'
    write(6,'(a)') '-dw          generate rmc distance windows file, requires file input'
    write(6,'(a)') '-elabel      generate new element based labels'
    write(6,'(a)') '-field       write bonds part of DLPOLY FIELD file, requires file input'
    write(6,'(a)') '-gasp        write the configuration in GASP input format'
    write(6,'(a)') '-gridcell    read size of grid cell to replace the default value of 5.0'
    write(6,'(a)') '-gulpc       write structure in format for GULP, in cartesian coordinates'
    write(6,'(a)') '-gulpf       write structure in format for GULP, in fractional coordinates'
    write(6,'(a)') '-hollow      creates hollow nanoparticles'
    write(6,'(a)') '-history     list version history'
    write(6,'(a)') '-keep        Keep orginal coordinates, ie so not shift atom coordinates to put into unit cell'
    write(6,'(a)') '-limits      Restrict calculations to a range of target atoms specified within [...] brackets'
    write(6,'(a)') '-listf       Provide a list file as an argument'
    write(6,'(a)') '-metadata    read metadata from a file containing keyword/value pairs'
    write(6,'(a)') '-molecules   Find molecules in a configuration, using file with provided name'
    write(6,'(a)') '-nano        create a nanoparticle with dimensions specified within [...] brackets as argument *'
    write(6,'(a)') '-nanoscale   Scale the size of a nanoparticle with factor specified within [...] brackets as argument'
    write(6,'(a)') '-nanobox     Specify the multipliers for the nanoparticle box. If not provided, [4 4 4] will be used.'
    write(6,'(a)') '-noannotate  request no annotations to files'
    write(6,'(a)') '-not         do not generate the Bragg scattering files if the .exp and .lst files are present'
    write(6,'(a)') '-one         only one unit cell required'
    write(6,'(a)') '-onecell     use only one unit cell for analysis'
    write(6,'(a)') '-order       order the atoms by user-supplied list in [...] brackets *'
    write(6,'(a)') '-orient      rotate the molecules randomly according to the keyword given as argument. Options are:'
    write(6,'(a)') '             random, orth or 90, flipx or 180x, flipy or 180y, flipz or 180z, flipr or 180r' 
    write(6,'(a)') '-origin      Reset of the configuration origin based on the setting the mean position at the centre'
    write(6,'(a)') '-output      provide an output directory/folder as supplied argument'
    write(6,'(a)') '-pdf         provide a calculated pdf with maximum distance as supplied argument'
    write(6,'(a)') '-pdfwidth    provide a width for broadening the pdf peaks as supplied argument'
    write(6,'(a)') '-reduce      move all atoms in a supercell back into the original unit cell'
    write(6,'(a)') '-relabel     generate new integer-based atom labels'
    write(6,'(a)') '-reorder     reorder atoms without generating supercell'
    write(6,'(a)') '-rmc3        write RMC v3 configuration file (extension .cfg)'
    write(6,'(a)') '-rmc6f       write RMC v6f configuration file'
    write(6,'(a)') '-rmc6o       write RMC v6o configuration file'
    write(6,'(a)') '-rect        create rectangle cell from hexagonal'
    write(6,'(a)') '-rmcanal     analyse rmc6f file based on the RMCprofile molecules files'
    write(6,'(a)') '-replace     replace atoms as provided within [...] brackets'
    write(6,'(a)') '-rsame       change the default test parameter to determine identical sites'
    write(6,'(a)') '-shape       provide nanoparticle shape as argument (e.g. sphere, cube, cylinder, rectangle, ellipsoid etc.)'
    write(6,'(a)') '-shift       makes sure all fractional coordinate values are between 0 and 1'
    write(6,'(a)') '-shiftx      shift all coordinates along x direction by the amount given as a number. Also y and z'
    write(6,'(a)') '-silica      convert a www file input from silicon into silica'
    write(6,'(a)') '-size        tune supercell to approximate given size, requires argument'
    write(6,'(a)') '-sort        sort the atoms by atom type'
    write(6,'(a)') '-supercell   provide the supercell multipliers within [...] brackets as argument *'
!    write(6,'(a)') '-transform   provide a transformation matrix within [...] brackets as argument'
    write(6,'(a)') '-Uiso        apply Gaussian distribution displacement for all atoms with Uiso value specified in [...]'
    write(6,'(a)') '-Biso        apply Gaussian distribution displacement for all atoms with Biso value spacified in [...]'
    write(6,'(a)') '-usehist     when reading a GSAS EXP file, use only the histograms specified within [...] brackets'
    write(6,'(a)') '-uselabels   use atom labels as provided in structure file'
    write(6,'(a)') '-vacancy     create vacancies, provide atoms and fractions within [...] brackets as argument'
    write(6,'(a)') '-version     prints the version number'
    write(6,'(a)') '-writef      in some cases where you have an option, write the output file in fractional coordinates'
    write(6,'(a)') '-writeo      in some cases where you have an option, write the output file in orthogonal coordinates'
    write(6,'(a)') '-www         write out a file with tag .www in same format used in the WWW program'
    write(6,'(a)') '-xtl         write an XTL format configuration file'
    write(6,'(a)') '-ylm         Analyse orientational variables, and nocentre (or nocenter) after to exclude central atom'
    write(6,'(a)') '* denotes that if the [...] brackets are not provided the information can be given interactvely'
    write(6,'(a)') '------------------------'
    write(6,'(a)') 'File types are defined by their extension; allowed extensions are:'
    write(6,'(a)') '.TBL or .tbl     - files produced by GSAS'
    write(6,'(a)') '.CIF or .cif     - crystal information format'
    write(6,'(a)') '.CELL or .cell   - Discus format'
    write(6,'(a)') '.SFF or .sff     - our own simple file format'
    write(6,'(a)') '.CFG or .cfg     - RMC v3 configuration file'
    write(6,'(a)') '.HIS or .his     - RMC v3 histogram file'
    write(6,'(a)') '.RES or .res     - old SHELX structure file, as used by CASTEP'
    write(6,'(a)') '.RMC6F or .rmc6f - RMC v6 configuration file'
    write(6,'(a)') '.CSSR or .cssr   - CSSR file'
    write(6,'(a)') '.WWW or .www     - generated by the Vink WWW program'
!    write(6,'(a)') '.CML or .cml     - Chemical Markup Language file *'
!    write(6,'(a)') '.XML or .xml     - assumed CML language *'
    write(6,'(a)')
    write(6,'(a)') ' DL_POLY''s CONFIG and REVCON files are supported without file name extension'
    write(6,'(a)') '------------------------'
!    write(6,'(a)') '* denotes this functionality is not yet implemented but will be soon'
    stop
    end if

    ncell = 1

    do i = 1,nargs
      call get_command_argument(i,temp)
      temp = adjustl(temp) ; arg(i) = temp(1:max_arg_length)

      if (temp(1:8)=='-version') then
         write(6,'(a)') 'data2config, version '//trim(version_number(nversions)) &
           //' ('//trim(version_date(nversions))//')'
        write(6,'(a/)') 'Last change is: '//trim(version_history(nversions))
        stop
      end if

      if (temp(1:8)=='-history') then
        write(6,'(a)') 'data2config, version history'
        write(6,'(a)') 'current version is: '//trim(version_number(nversions))
        write(6,'(a)') 'last update is: '//trim(version_history(nversions))
        do n = 1,nversions
          write(6,'(a,a,a)') version_number(n),version_date(n),': '//trim(version_history(n))
        end do
        stop
      end if

      if (temp(1:5)=='-help') call help_info

      if (temp(1:5)=='-cssr') lcssr = .true.

      if (temp(1:5)=='-cmtx') lcmtx = .true.

      if (temp(1:4)=='-cif') lcif = .true.

      if (temp(1:7)=='-castep') lcastep = .true.

      if (temp(1:7)=='-dlpoly') ldlpoly = .true.
      
      if (temp(1:4)=='-www') lwww_write = .true.

      if (temp(1:6)=='-gulpf') then
        lgulp = .true.
      end if

      if (temp(1:7)=='-writef') then
        lwritef = .true.
        lwriteo = .false.
      end if

      if (temp(1:7)=='-writeo') then
        lwritef = .false.
        lwriteo = .true.
      end if

      if (temp(1:6)=='-gulpc') then
        lgulp = .true.
        lgulpc = .true.
      end if

      if (temp(1:5)=='-keep') ladjust_xyz_signs = .false.

      if (temp(1:6)=='-field') then
        lfield = .true.
        call get_command_argument(i+1,fieldfile)
        fieldfile = adjustl(fieldfile)
        inquire(file=trim(fieldfile),exist=lokay)
        write(6,'(a)') 'data2config, version '//trim(version_number(nversions))
        write(6,'(a)') 'Fieldfile = '//trim(fieldfile)
        if (.not.lokay) then
          write(6,'(a)') 'Fieldfile not found'
          stop
        end if
      end if

      if (temp(1:10)=='-molecules') then
        lmolecules = .true.
        call get_command_argument(i+1,moleculefile)
        moleculefile = adjustl(moleculefile)
        inquire(file=trim(moleculefile),exist=lokay)
        write(6,'(a)') 'data2config, version '//trim(version_number(nversions))
        write(6,'(a)') 'Moleculefile = '//trim(moleculefile)
        if (.not.lokay) then
          write(6,'(a)') 'Moleculefile not found'
          stop
        end if
      end if

      if (temp(1:3)=='-dw') then
        ldw = .true.
        call get_command_argument(i+1,dwfile)
        dwfile = adjustl(dwfile)
        inquire(file=trim(dwfile),exist=lokay)
        write(6,'(a)') 'data2config, version '//trim(version_number(nversions))
        write(6,'(a)') 'Distance Windows file = '//trim(dwfile)
        if (.not.lokay) then
          write(6,'(a)') 'Distance Windows file not found'
          stop
        end if
      end if

      if (temp(1:11)=='-potentials') then
        lpotentials = .true.
        call get_command_argument(i+1,potentialsfile)
        potentialsfile = adjustl(potentialsfile)
        inquire(file=trim(potentialsfile),exist=lokay)
        write(6,'(a)') 'data2config, version '//trim(version_number(nversions))
        write(6,'(a)') 'Potentialsfile = '//trim(potentialsfile)
        if (.not.lokay) then
          write(6,'(a)') 'Potentialsfile not found'
          stop
        end if
      end if

      if ((temp(1:5)=='-size').and..not.lnano) then
        call get_command_argument(i+1,clength)
        ierror = 0
        read(clength,*,iostat=ierror) approxsize
        if (ierror/=0) then
          write(6,'(a)') 'data2config, version '//trim(version_number(nversions))
          write(6,'(a)') 'Error on input line; expected argument giving approximate cell size'
          stop
        end if
        lsize = .true.
      end if

      if (temp(1:7)=='-shiftx') then
        call get_command_argument(i+1,clength)
        ierror = 0
        read(clength,*,iostat=ierror) shiftx
        if (ierror/=0) then
          write(6,'(a)') 'data2config, version '//trim(version_number(nversions))
          write(6,'(a)') 'Error on input line; expected argument giving the size of a shift'
          stop
        end if
        lshiftx = .true.
      end if

      if (temp(1:7)=='-shifty') then
        call get_command_argument(i+1,clength)
        ierror = 0
        read(clength,*,iostat=ierror) shifty
        if (ierror/=0) then
          write(6,'(a)') 'data2config, version '//trim(version_number(nversions))
          write(6,'(a)') 'Error on input line; expected argument giving the size of a shift'
          stop
        end if
        lshifty = .true.
      end if

      if (temp(1:7)=='-shiftz') then
        call get_command_argument(i+1,clength)
        ierror = 0
        read(clength,*,iostat=ierror) shiftz
        if (ierror/=0) then
          write(6,'(a)') 'data2config, version '//trim(version_number(nversions))
          write(6,'(a)') 'Error on input line; expected argument giving the size of a shift'
          stop
        end if
        lshiftz = .true.
      end if

      if (temp(1:7)=='-origin') lorigin = .true.

      if (temp(1:9)=='-gridcell') then
        call get_command_argument(i+1,clength)
        ierror = 0
        read(clength,*,iostat=ierror) gridcellsize
        if (ierror/=0) then
          write(6,'(a)') 'data2config, version '//trim(version_number(nversions))
          write(6,'(a)') 'Error on input line; expected argument giving grid cell size'
          stop
        end if
      end if

      if (temp(1:8)=='-compare') then
        lcompare = .true.
        call get_command_argument(i+1,comparefile)
        comparefile = adjustl(comparefile)
        inquire(file=trim(comparefile),exist=lokay)
        write(6,'(a)') 'data2config, version '//trim(version_number(nversions))
        write(6,'(a)') 'File for comparison = '//trim(comparefile)
        if (.not.lokay) then
          write(6,'(a)') 'Comparison file not found'
          stop
        end if
      end if

      if (temp(1:6)=='-align') then
        lalign = .true.
        call get_command_argument(i+1,alignfile)
        alignfile = adjustl(alignfile)
        inquire(file=trim(alignfile),exist=lokay)
        write(6,'(a)') 'data2config, version '//trim(version_number(nversions))
        write(6,'(a)') 'Reference file for alignment = '//trim(alignfile)
        if (.not.lokay) then
          write(6,'(a)') 'Reference file for alignment not found'
          stop
        end if
      end if

      if (temp(1:7)=='-reduce') lreduce = .true.

      if (temp(1:4)=='-pdf'.and.(.not.lpdf)) then
        call get_command_argument(i+1,clength)
        ierror = 0
        read(clength,*,iostat=ierror) rpdfmax
        if (ierror/=0) then
          write(6,'(a)') 'data2config, version '//trim(version_number(nversions))
          write(6,'(a)') 'Error on input line; expected argument giving pdf size limit'
          stop
        end if
        lpdf = .true.
      end if

      if (temp(1:9)=='-pdfwidth') then
        call get_command_argument(i+1,clength)
        ierror = 0
        read(clength,*,iostat=ierror) pdfwidth
        if (ierror/=0) then
          write(6,'(a)') 'data2config, version '//trim(version_number(nversions))
          write(6,'(a)') 'Error on input line; expected argument giving pdf peak width'
          stop
        end if
        lpdfbroaden = .true.
      end if

      if (temp(1:6)=='-listf') then
        llistfile = .true.
        call get_command_argument(i+1,listfile)
        listfile = adjustl(listfile)
        inquire(file=trim(listfile),exist=lokay)
        write(6,'(a)') 'data2config, version '//trim(version_number(nversions))
        write(6,'(a)') 'PDF origins list file = '//trim(listfile)
        if (.not.lokay) then
          write(6,'(a)') 'PDF origins list file not found'
          stop
        end if
      end if

! Note that the -alpha input is planned to be an undocumented option, allowing passing
! of a single variable for any customised versions of data2config
      if (temp(1:5)=='-zeta') then
        call get_command_argument(i+1,clength)
        ierror = 0
        read(clength,*,iostat=ierror) zeta
        if (ierror/=0) then
          write(6,'(a)') 'data2config, version '//trim(version_number(nversions))
          write(6,'(a)') 'Error on input line; expected argument giving value of zeta'
          stop
        end if
        lzeta = .true.
      end if
      
      if (temp(1:8)=='-relabel') lrelabel = .true.

      if (temp(1:8)=='-elabel') lelabel = .true.

      if (temp(1:10)=='-uselabels') luselabels = .true.
      
      if ((temp(1:10)=='-supercell').and..not.lnano) then
        lsupercell = .true.
        n = index(comline,'-supercell') + 10 
        ctemp = adjustl(comline(n:))
        n1 = index(ctemp,'[') ; n2 = index(ctemp(n1+1:),']')
        if ((n1>0).and.(n2>0)) then
          lsupercell_list = .true.
          ierror = 0
          read(ctemp(n1+1:n2),*,iostat=ierror) ncell
          if (ierror/=0) lsupercell_list = .false.
        end if
      end if

      if (temp(1:7)=='-limits') then
        lxyzlimits = .true.
        n = index(comline,'-limits') + 7 
        ctemp = adjustl(comline(n:))
        n1 = index(ctemp,'[') ; n2 = index(ctemp(n1+1:),']')
        if ((n1>0).and.(n2>0)) then
          ierror = 0
          read(ctemp(n1+1:n2),*,iostat=ierror) xyzlimits
          if (ierror/=0) then
            write(6,'(a)') '6 values to define the limits of the x, y and z fractional coordinates were expected'
            stop
          end if
          if (xyzlimits(1)>=xyzlimits(2)) then
            write(6,'(a)') 'xmax should be greater than xmin, but it seems not to be'
            stop
          end if
          if (xyzlimits(3)>=xyzlimits(4)) then
            write(6,'(a)') 'ymax should be greater than ymin, but it seems not to be'
            stop
          end if
          if (xyzlimits(5)>=xyzlimits(6)) then
            write(6,'(a)') 'zmax should be greater than zmin, but it seems not to be'
            stop
          end if
        end if
      end if

      if ((temp(1:10)=='-transform').and..not.lnano) then
        n = index(comline,'-transform') + 10 
        ctemp = adjustl(comline(n:))
        n1 = index(ctemp,'[') ; n2 = index(ctemp(n1+1:),']')
        if ((n1==0).or.(n2==0)) then
          write(6,'(a)') 'data2config, version '//trim(version_number(nversions))
          write(6,'(a)') 'We expected a 3x3 transformation matrix written as nine integers within'
          write(6,'(a)') '[...] brackets. These were not found, so we assume an error'
          stop
        end if
        read(ctemp(n1+1:n2),*) mtransform
        ltransform = .true.
      end if

      if (temp == '-usehist') then
        n = index(comline,'-usehist') + 8
        ctemp = adjustl(comline(n:))
        n1 = index(ctemp,'[') ; n2 = index(ctemp(n1+1:),']')
        if ((n1>0).and.(n2>0)) then
          ierror = 0
          nusehist = 0
          ctempcount = ctemp(n1+1:n1+n2-1)

          read(ctempcount, *, iostat=ierror) iusehist
          do while (ierror == 0)
             nusehist = nusehist + 1
             ctempcount = adjustl(ctempcount(index(ctempcount, ' '):))
             read(ctempcount, *, iostat=ierror) iusehist
          end do

          allocate(useonly(nusehist))
          read(ctemp(n1+1:n2), *) useonly
        end if
      end if

      if (temp(1:8)=='-vacancy') then
        n = index(comline,'-vacancy') + 8 
        ctemp = adjustl(comline(n:))
        n1 = index(ctemp,'[') ; n2 = index(ctemp(n1+1:),']')
        if ((n1==0).or.(n2==0)) then
          write(6,'(a)') 'data2config, version '//trim(version_number(nversions))
          write(6,'(a)') 'We expected information about vacancy atoms given within'
          write(6,'(a)') '[...] brackets. These were not found, so we assume an error'
          stop
        end if
        cvacancy_list = trim(adjustl(ctemp(n1+1:n2)))
        lvacancy = .true.
      end if

      if (temp(1:7)=='-delete') ldelete = .true.

      if (temp(1:8)=='-replace') then
        n = index(comline,'-replace') + 8 
        ctemp = adjustl(comline(n:))
        n1 = index(ctemp,'[') ; n2 = index(ctemp(n1+1:),']')
        if ((n1==0).or.(n2==0)) then
          write(6,'(a)') 'data2config, version '//trim(version_number(nversions))
          write(6,'(a)') 'We expected information about replacement atoms given within'
          write(6,'(a)') '[...] brackets. These were not found, so we assume an error'
          stop
        end if
        creplacement_list = trim(adjustl(ctemp(n1+1:n2)))
        lreplace = .true.
      end if

      if (temp(1:6)=='-rsame') then
        call get_command_argument(i+1,clength)
        ierror = 0
        read(clength,*,iostat=ierror) rsame
        if (ierror/=0) then
          write(6,'(a)') 'data2config, version '//trim(version_number(nversions))
          write(6,'(a)') 'Error on input line; expected argument giving size test parameter'
          stop
        end if
      end if

      rotate: if (temp(1:7)=='-rotate') then
       lrotate = .true.
       n = index(comline,'-rotate') + 7 
        ctemp = adjustl(comline(n:)) !read all line and assign "value" to ctemp
        n1 = index(ctemp,'[') ; n2 = index(ctemp(n1+1:),']') !assign value '[' to n1 and n2 is ']'
        if ((n1>0).and.(n2>0)) then !n1>0 if [ existes in the string otherwise is zero
          lrotate_list = .true.
          ierror = 0
          read(ctemp(n1+1:n2),*,iostat=ierror) fi,csi,psi 
           
          if ((ierror/=0)) lrotate_list = .false.
         end if
      end if rotate
       
      nano: if (temp(1:5)=='-nano') then
        lnano = .true.
        n = index(comline,'-nano') + 5 
        ctemp = adjustl(comline(n:)) !read all line and assign "value" to ctemp
        n1 = index(ctemp,'[') ; n2 = index(ctemp(n1+1:),']') !assign value '[' to n1 and n2 is ']'
        if ((n1>0).and.(n2>0)) then !n1>0 if [ existes in the string otherwise is zero
          lnano_list = .true.
          ierror = 0
          read(ctemp(n1+1:n2),*,iostat=ierror) rnano !read numbers in between n1+1 and n2 and assign the values to rnano
         if ((ierror/=0).or.(rnano(1)<=0.or.rnano(2)<=0.or.rnano(3)<=0)) lnano_list = .false.
        else
          read(ctemp,*,iostat=ierror) rnano(1)
          if (ierror==0) then
            rnano = rnano(1)
            lnano_list = .true.
          else
            lnano_list = .false.
          end if
        end if
        lsize = .true.
        approxsize =  max(rnano(1),rnano(2),rnano(3))

        ! >>>>>>>>>>>> Yuanpeng -> Add in the functionality for user to supply >>>>>>>>
        ! the box size for the nanoparticle. >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        n = index(comline,'-nanobox') + 8
        ctemp = adjustl(comline(n:))
        n1 = index(ctemp,'[') ; n2 = index(ctemp(n1+1:),']')
        if ((n1>0).and.(n2>0)) then
          lnanobox_list = .true.
          ierror = 0
          read(ctemp(n1+1:n2),*,iostat=ierror) nanobox
         if ((ierror/=0).or.(nanobox(1)<=0.or.nanobox(2)<=0.or.nanobox(3)<=0)) lnanobox_list = .false.
        end if
        ! <<<<<<<<<<<< Yuanpeng finishes here <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<             
      end if nano

      if (temp(1:10)=='-nanoscale') then
        call get_command_argument(i+1,clength)
        ierror = 0
        read(clength,*,iostat=ierror) scalenano
        if (ierror/=0) then
          write(6,'(a)') 'data2config, version '//trim(version_number(nversions))
          write(6,'(a)') 'Error on input line; expected argument giving value of scalenano'
          stop
        end if
        lnanoscale = .true.
      end if

      orient: if (temp(1:7)=='-orient') then
        call get_command_argument(i+1,clength)
        lorient = .true.
        ierror = 0
        read(clength,*,iostat=ierror) oshape
        select case(trim(oshape))
          case('orth','90','random','flipx','flipy','flipz','flipr','180x','180y','180z','180r')
            continue
          case default
            write(6,'(a)') 'Orient shape '//trim(oshape)//' not supported'
            stop
         end select 
        if (ierror/=0) then
          write(6,'(a)') 'data2config, version '//trim(version_number(nversions))
          write(6,'(a)') 'Error on input line; expected argument giving orient shape'
          stop
        end if
      end if orient
        
     hollow: if (temp(1:7)=='-hollow') then
     lhollow = .true.
         n = index(comline,'-hollow') + 7 
        ctemp = adjustl(comline(n:)) !read all line and assign "value" to ctemp
        n1 = index(ctemp,'[') ; n2 = index(ctemp(n1+1:),']') !assign value '[' to n1 and n2 is ']'
        if ((n1>0).and.(n2>0)) then !n1>0 if [ existes in the string otherwise is zero
      lhollow_list = .true.
          ierror = 0
          read(ctemp(n1+1:n2),*,iostat=ierror) rhollow !read numbers in between n1+1 and n2 and assign the values to rhollow
        
          if ((ierror/=0).or.(rhollow(1)<=0.or.rhollow(2)<=0.or.rhollow(3)<=0)) lhollow_list = .false.
         
      end if
             
     end if hollow
          
      
      if (temp(1:6)=='-shape') then
        call get_command_argument(i+1,clength)
        cshape = trim(clength)
        select case(trim(cshape))
          case('sphere')
            continue
          case('cube')
            continue
          case('cylinder')
            continue
          case('rectangle')
            continue
          case('ellipsoid')
            continue
          case default
            write(6,'(a)') 'Nanoparticle shape '//trim(cshape)//' not supported'
            stop 
        end select
      end if

      if (temp(1:9)=='-metadata') then
        call get_command_argument(i+1,fmetadata)
        lmetadata = .true. ; lannotate = .false.
        inquire(file=trim(fmetadata),exist=lokay)
        write(6,'(a)') 'data2config, version '//trim(version_number(nversions))
        if (.not.lokay) then
          write(6,'(a)') 'Metadata file not found'
          stop
        end if
      end if 

      if (temp(1:8)=='-cmbonds') then
        call get_command_argument(i+1,cmbondfile)
        lcrystalmakerbonds = .true.
        inquire(file=trim(cmbondfile),exist=lokay)
        if (.not.lokay) then
          write(6,'(a)') 'CrystalMaker bonds file not found'
          stop
        end if
      end if 

      if (temp(1:7)=='-cmsort') then
        call get_command_argument(i+1,cmbondfile)
        inquire(file=trim(cmbondfile),exist=lokay)
        if (.not.lokay) then
          write(6,'(a)') 'CrystalMaker bonds file not found'
          stop
        end if
        lcmsort = .true.
      end if 

      if (temp(1:7)=='-output') then
        call get_command_argument(i+1,wfolder)
        n = len_trim(wfolder)
        if (lunix.and.(wfolder(n:n)/='/')) wfolder = trim(wfolder)//'/'
        if (lwindows.and.(wfolder(n:n)/='\')) wfolder = trim(wfolder)//'\'
        inquire(file=trim(wfolder),exist=lokay)
        write(6,'(a)') 'data2config, version '//trim(version_number(nversions))
        if (.not.lokay) then
          write(6,'(a)') 'Output folder not found'
          stop
        end if
      end if

      if (temp(1:5)=='-bank') then
        call get_command_argument(i+1,clength)
        ierror = 0
        read(clength,*,iostat=ierror) mbank
        if (ierror/=0) then
          write(6,'(a)') 'data2config, version '//trim(version_number(nversions))
          write(6,'(a)') 'Error on input line; expected argument giving bank number'
          stop
        end if
      end if


      if (temp(1:4)=='-cml') lcml = .true.

      if (temp(1:5)=='-rmc3') lrmc3 = .true.

      if (temp(1:6)=='-rmc6f') lrmc6f = .true.

      if (temp(1:6)=='-rmc6o') lrmc6o = .true.

      if (temp(1:8)=='-crystal') lcrystal = .true.

      if (temp(1:6)=='-crush') lcrush = .true.

      if (temp(1:5)=='-rect') lrect = .true.

      if (temp(1:5)=='-sort') lsort = .true.
      
      if (temp(1:8)=='-atomeye') latomeye = .true.

      if (temp(1:5)=='-gasp') lgasp = .true.

      if (temp(1:7)=='-silica') lsilica = .true.

      if (temp(1:4)=='-xtl') lxtl = .true.
      
      if (temp(1:5)=='-Uiso') then
        lUiso = .true.
        n = index(comline,'-Uiso') + 5
        ctemp = adjustl(comline(n:))
        n1 = index(ctemp,'[') ; n2 = index(ctemp(n1+1:),']')
        n3 = index(ctemp,'-')
        if (n3 == 0) n3 = len(ctemp)
        if ((n1>0).and.(n2>0).and.(n3>n1)) then
          cUiso_list = trim(adjustl(ctemp(n1+1:n2)))
          lUiso_list = .true.
          Uiso_factor = 1.0
        end if
      end if   
      
      if (temp(1:5)=='-Biso') then
        lUiso = .true.
        n = index(comline,'-Biso') + 5
        ctemp = adjustl(comline(n:))
        n1 = index(ctemp,'[') ; n2 = index(ctemp(n1+1:),']')
        n3 = index(ctemp,'-')
        if (n3 == 0) n3 = len(ctemp)
        if ((n1>0).and.(n2>0).and.(n3>n1)) then
          cUiso_list = trim(adjustl(ctemp(n1+1:n2)))
          lUiso_list = .true.
          Uiso_factor = 1.0 / (8.0 * pi * pi)
        end if
      end if 

      if (temp(1:6)=='-order') then
        lorder = .true.
        n = index(comline,'-order') + 6 
        ctemp = adjustl(comline(n:))
        n1 = index(ctemp,'[') ; n2 = index(ctemp(n1+1:),']')
        if ((n1>0).and.(n2>0)) then
          corder_list = trim(adjustl(ctemp(n1+1:n2)))
          lorder_list = .true.
        end if
      end if

      if (temp(1:8)=='-reorder') then
        lorder = .true.
        n = index(comline,'-reorder') + 8 
        ctemp = adjustl(comline(n:))
        n1 = index(ctemp,'[') ; n2 = index(ctemp(n1+1:),']')
        if ((n1>0).and.(n2>0)) then
          corder_list = trim(adjustl(ctemp(n1+1:n2)))
          lorder_list = .true.
          lreorder = .true.
          lannotate = .false.
          lone = .true.
        end if
      end if

      if (temp(1:5)=='-diag') ldiag = .true.

      if (temp(1:4)=='-not') lnot = .true.  ! Undocumented fudge

      if (temp(1:4)=='-one') lone = .true.

      if (temp(1:8)=='-onecell') lonecell = .true.

      if (temp(1:9)=='-analysis') then
        lanal = .true.
        call get_command_argument(i+1,bondfile)
        bondfile = adjustl(bondfile)
        write(6,'(a)') 'Bondfile = '//trim(bondfile)
      end if
      
      if (temp(1:7)=='-angles') langles = .true.

      if (temp(1:4)=='-ylm') then
        lylm = .true.
        call get_command_argument(i+1,buffer)
        buffer = adjustl(buffer)
        lmolcentre = .true.
        if (buffer(1:3)=='noc')  lmolcentre = .false.
      end if
      
      if (temp(1:11)=='-noannotate') lannotate = .false.

      if (temp(1:4)=='-mag') then
        lmag = .true.
        call get_command_argument(i+1,filemag)
        filemag = adjustl(filemag)
      end if

      if (temp(1:7)=='-dlanal') ldlanal = .true.

      if (temp(1:8)=='-rmcanal') lrmcanal = .true.

    end do

    lfield = (lfield.and.ldlpoly)
    call get_command_argument(nargs,filename)
    filename = adjustl(filename)
    if (index(filename,'\')>0) then
      lwindows = .true. ; lunix = .false.
    end if
    fileroot = trim(filename)
    n = index(filename,'\',back=.true.)
    if (n>0) fileroot = filename(n+1:)
    n = index(filename,'/',back=.true.)
    if (n>0) fileroot = filename(n+1:)

!    if (lrmc3) lsort = .true.
!    if (lcrystal) lsort = .true.
    if (lcrystal) lannotate = .false.
    if (lcrush) lannotate = .false.
    if (lrmc6f.and.(.not.lorder).and.(.not.lsort).and.(.not.lalign)) lsort = .true.
    if (lrmc6o.and.(.not.lorder).and.(.not.lsort).and.(.not.lalign)) lsort = .true.

    if (trim(filename)=='no filename given') then
      write(6,'(a)') 'data2config, version '//trim(version_number(nversions))
      write(6,'(a)') char(9)//&
        'You need to give the file name as the parameter to this program'
      stop
    end if
    
   inquire(file='molecules.dat', exist = lexist)
   if (lylm.and.(.not.lexist).and.(.not.lmolecules)) then
      lylm = .false.
      write(6,'(a)') 'Request for orientation analysis ignored because no -molecules flag included'
    end if

! Set a flag for offering the option of a supercell
    if (lreduce) lone = .true.
    if (lcml) lsupercell = .true.
    if (lrmc3) lsupercell = .true.
    if (lrmc6f) lsupercell = .true.
    if (lrmc6o) lsupercell = .true.
    if (lcssr) lsupercell = .true.
    if (lcif) lsupercell = .true.
    if (ldlpoly) lsupercell = .true.
    if (lgulp) lsupercell = .true.
    if (lcrystal) lsupercell = .true.
    if (ltransform) lsupercell = .false.
    !if (lnano) lsupercell = .false. !ask
    if (lone) lsupercell = .false.
    if (lrmc6f) lrelabel = .true.
    if (lrmc6o) lrelabel = .true.
    if (lrelabel.and.luselabels) luselabels = .false.

! Here look for input file type by filename extension
    inquire(file=trim(filename),exist=lexist)
    if (.not.lexist) then
      write(6,'(a)') 'Main structure file does not appear to exist'
      stop
    end if
    n = index(filename,'.')
    lokay = .false.

    if ((trim(filename(n+1:))=='cml').or.(trim(filename(n+1:))=='CML') ) then
      cinput = 'cml'
      lokay = .true.
      lcml = .false.
    end if

    if ((trim(filename(n+1:))=='xml').or.(trim(filename(n+1:))=='XML') ) then
      cinput = 'cml'
      lokay = .true.
      lcml = .false.
    end if

    if ((trim(filename(n+1:))=='tbl').or.(trim(filename(n+1:))=='TBL') ) then
      cinput = 'tbl'
      lokay = .true.
      filetest = filename(1:n-1)//'.EXP'
      inquire(file=trim(filetest),exist=l1)
      filetest = filename(1:n-1)//'.exp'
      inquire(file=trim(filetest),exist=l2)
      lexp = (l1.or.l2)
    end if

    if ((trim(filename(n+1:))=='cif').or.(trim(filename(n+1:))=='CIF') ) then
      cinput = 'cif'
      lokay = .true.
    end if

    if ((trim(filename(n+1:))=='his').or.(trim(filename(n+1:))=='HIS') ) then
      cinput = 'his'
      lokay = .true.
      lrmc3 = .false.
      lhis = .true.
    end if

    if ((trim(filename(n+1:))=='rmc6f').or.(trim(filename(n+1:))=='RMC6F') ) then
      cinput = 'rmc6f'
      lokay = .true.
      lconfigin = .true.
      luselabels = .true.
      lrelabel = .false.
    lelabel = .false.
    end if

    if ((trim(filename(n+1:))=='cssr').or.(trim(filename(n+1:))=='CSSR') ) then
      cinput = 'cssr'
      lokay = .true.
      lcssr = .false.
    end if

    if ((trim(filename(n+1:))=='cssr').or.(trim(filename(n+1:))=='CSSR') ) then
      cinput = 'cssr'
      lokay = .true.
      lcssr = .false.
    end if

    if ((trim(filename(n+1:))=='cfg').or.(trim(filename(n+1:))=='CFG') ) then
!     Check which type of .cfg file we have; RMCprofile or AtomEye. The only way is to read
!     the first line
      lokay = .true.
      lconfigin = .true.
    end if

    if ((trim(filename(n+1:))=='sff').or.(trim(filename(n+1:))=='SFF') ) then
      cinput = 'sff'
      lokay = .true.
    end if

    if ((trim(filename(n+1:))=='www').or.(trim(filename(n+1:))=='WWW') ) then
      lwww = .true.
      cinput = 'www'
      lokay = .true.
    end if

    if ((trim(filename(n+1:))=='res').or.(trim(filename(n+1:))=='RES') ) then
      cinput = 'res'
      lokay = .true.
    end if

    if ((trim(filename(n+1:))=='cfgcom').or.(trim(filename(n+1:))=='CFGCOM') ) then
      cinput = 'crystal'
      lcrystal = .false.
      lokay = .true.
    end if

    if ((trim(filename(n+1:))=='cell').or.(trim(filename(n+1:))=='CELL') ) then
    cinput = 'cell'
    lconfigin = .true.
    lokay = .true.
    end if

    if (filename(1:6)=='CONFIG'.and..not.lokay) then
      cinput = trim(filename)
      lconfigin = .true.
      lokay = .true.
    end if

    if (filename(1:6)=='REVCON'.and..not.lokay) then
      cinput = trim(filename)
      lconfigin = .true.
      lokay = .true.
    end if

    if (filename(1:7)=='HISTORY'.and..not.lokay) then
      cinput = trim(filename)
      lconfigin = .true.
      lokay = .true.
    end if

    if (.not.lokay) then
      write(6,'(a)') char(9)//&
      'Input filetype not supported; it needs to be one of cml, tbl, cif, cfg, res, rmc6f, CONFIG, REVCON, HISTORY or cssr'
    stop
    end if
    
    if (latomeye.and.lrmc6f) then
      write(6,'(a)') &
'The rmc3f and atomeye configurations both end with .cfg, which will cause problems.'
      write(6,'(a)') &
'Please consider running the program twice, once for each output type.'
      stop
    end if
    
    if (lconfigin.and.lpdf) gridcellsize = rpdfmax

    return

  end subroutine get_arguments

!===============================================================================

  subroutine get_symmetry_tbl
! ===========================

!--------------------------------------------------------------------------
!
! Read the symmetry information from the TBL file
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use channel_numbers

    implicit none

    integer :: i,j,n,m,ml
    character(len=132) :: buffer

    buffer = ''
    do while (index(buffer,'Space group')==0)
    read(ic,'(a)') buffer
    end do

    read(ic,'(a)') buffer  ! Here we read the line that gives lattice info
    n = index(buffer,'centric')
    lcentric = (buffer(n-1:n-1)/='a')  ! true if the we have centric rather than acentric
    lattice_info(2) = 0 ; if (lcentric) lattice_info(2) = 1
    n = index(buffer,'-centered')
    centre = 'P' ; if (n>0) centre = buffer(n-1:n-1)
    system = 'Nothing'
    if(index(buffer,'cubic')>0) system = 'cubic'
    if(index(buffer,'tetragonal')>0) system = 'tetragonal'
    if(index(buffer,'orthorhombic')>0) system = 'orthorhombic'
    if(index(buffer,'trigonal')>0) system = 'trigonal'
    if(index(buffer,'hexagonal')>0) system = 'hexagonal'
    if(index(buffer,'rhombohedral')>0) system = 'rhombohedral'
    if(index(buffer,'monoclinic')>0) system = 'monoclinic'
    if(index(buffer,'triclinic')>0) system = 'triclinic'
    read(ic,'(a)') buffer  ! Here we read the line that gives multiplicity
    read(buffer(36:),*) multiplicity
    nsym = multiplicity
    if (lcentric) nsym = nsym/2
    if (centre=='I') nsym = nsym/2
    if (centre=='A') nsym = nsym/2
    if (centre=='B') nsym = nsym/2
    if (centre=='C') nsym = nsym/2
    if (centre=='F') nsym = nsym/4
    if (centre=='R') nsym = nsym/3
    allocate(coperators(nsym))

    buffer = ''
    do while (index(buffer,'The equivalent positions are:')==0)
    read(ic,'(a)') buffer
    end do
    read(ic,*)
    ml = 3 ; n = nsym/ml ; m = modulo(nsym,ml)
    if (m>0) then
     ml = 2 ; n = nsym/ml ; m = modulo(nsym,ml)
    end if
    if (m>0) then
    ml = 1 ; n = nsym
    end if
    do i = 1,n
    read(ic,'(a)') buffer
    if (i<n) then
      read(buffer,'(3a26)') (coperators(ml*(i-1)+j),j=1,ml)
    else
      read(buffer,'(3a26)') (coperators(j),j=ml*(n-1)+1,nsym)
    end if
    end do

    lattice_info(1) = nsym
    do i = 1,nsym ; coperators(i) = adjustl(coperators(i)) ; end do
    do i = 1,nsym
    if (coperators(i)(1:1)=='(') then
      n = index(coperators(i),')')
      coperators(i) = coperators(i)(n+1:)
    end if
    end do

    if (ldiag) then
    write(main_output,'(a)') 'Output from get_symmetry_tbl'
    write(main_output,'(a/)') '============================'
    do i = 1,nsym
      write(main_output,'(a)') trim(coperators(i))
    end do
    write(main_output,'(/)')
    end if

    return

  end subroutine get_symmetry_tbl

!===============================================================================

  subroutine get_structure_tbl
! ============================

!--------------------------------------------------------------------------
!
! Reads the structural details from the TBL file
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use channel_numbers

    implicit none

    integer :: i,j,n1,n2,ios,n
    character(len=132) :: buffer
    character(len=max_atom_string), allocatable :: caxyzt(:)
    double precision, allocatable :: occt(:),xft(:),yft(:),zft(:)
    character(len=24), allocatable :: citem(:)
    character(len=24), allocatable :: citem_temp(:)
    character(len=4), allocatable :: atom_namet(:)

! Now extract the lattice parameters
    ios = 0
    buffer = ''
    do while (index(buffer,'Lattice constants are')==0)
      read(ic,'(a)') buffer
    end do
    read(ic,'(a)') buffer ! Read lattice parameters
    if (trim(system)=='cubic') then
      n1 = index(buffer,'=') + 1 ; buffer = buffer(n1:) ; n2 = index(buffer,'(') - 1
      read(buffer(1:n2),*) a ; b = a ; c = a
      alpha = 90.0d0 ; beta = alpha ; gamma = alpha
    else if (trim(system)=='tetragonal') then
      n1 = index(buffer,'=') + 1 ; buffer = buffer(n1:) ; n2 = index(buffer,'(') - 1
      read(buffer(1:n2),*) a ; b = a
      n1 = index(buffer,'=') + 1 ; buffer = buffer(n1:)
      n1 = index(buffer,'=') + 1 ; buffer = buffer(n1:) ; n2 = index(buffer,'(') - 1
      read(buffer(1:n2),*) c
      alpha = 90.0d0 ; beta = alpha ; gamma = alpha
    else if (trim(system)=='orthorhombic') then
      n1 = index(buffer,'=') + 1 ; buffer = buffer(n1:) ; n2 = index(buffer,'(') - 1
      read(buffer(1:n2),*) a
      n1 = index(buffer,'=') + 1 ; buffer = buffer(n1:) ; n2 = index(buffer,'(') - 1
      read(buffer(1:n2),*) b
      n1 = index(buffer,'=') + 1 ; buffer = buffer(n1:) ; n2 = index(buffer,'(') - 1
      read(buffer(1:n2),*) c
      alpha = 90.0d0 ; beta = alpha ; gamma = alpha
    else if (trim(system)=='hexagonal') then
      n1 = index(buffer,'=') + 1 ; buffer = buffer(n1:) ; n2 = index(buffer,'(') - 1
      read(buffer(1:n2),*) a ; b = a
      n1 = index(buffer,'=') + 1 ; buffer = buffer(n1:)
      n1 = index(buffer,'=') + 1 ; buffer = buffer(n1:) ; n2 = index(buffer,'(') - 1
      read(buffer(1:n2),*) c
      alpha = 90.0d0 ; beta = alpha ; gamma = 120.0d0
    else if (trim(system)=='trigonal') then
      n1 = index(buffer,'=') + 1 ; buffer = buffer(n1:) ; n2 = index(buffer,'(') - 1
      read(buffer(1:n2),*) a ; b = a
      n1 = index(buffer,'=') + 1 ; buffer = buffer(n1:)
      n1 = index(buffer,'=') + 1 ; buffer = buffer(n1:) ; n2 = index(buffer,'(') - 1
      read(buffer(1:n2),*) c
      alpha = 90.0d0 ; beta = alpha ; gamma = 120.0d0
    else if (trim(system)=='monoclinic') then
      n1 = index(buffer,'=') + 1 ; buffer = buffer(n1:) ; n2 = index(buffer,'(') - 1
      read(buffer(1:n2),*) a
      n1 = index(buffer,'=') + 1 ; buffer = buffer(n1:) ; n2 = index(buffer,'(') - 1
      read(buffer(1:n2),*) b
      n1 = index(buffer,'=') + 1 ; buffer = buffer(n1:) ; n2 = index(buffer,'(') - 1
      read(buffer(1:n2),*) c
      read(ic,'(a)') buffer
      n1 = index(buffer,'=') + 1 ; buffer = buffer(n1:)
      n1 = index(buffer,'=') + 1 ; buffer = buffer(n1:) ; n2 = index(buffer,'(') - 1
      read(buffer(1:n2),*) beta
      alpha = 90.0d0 ; gamma = 90.0d0
    else if (trim(system)=='triclinic') then
      n1 = index(buffer,'=') + 1 ; buffer = buffer(n1:) ; n2 = index(buffer,'(') - 1
      read(buffer(1:n2),*) a
      n1 = index(buffer,'=') + 1 ; buffer = buffer(n1:) ; n2 = index(buffer,'(') - 1
      read(buffer(1:n2),*) b
      n1 = index(buffer,'=') + 1 ; buffer = buffer(n1:) ; n2 = index(buffer,'(') - 1
      read(buffer(1:n2),*) c
      read(ic,'(a)') buffer
      n1 = index(buffer,'=') + 1 ; buffer = buffer(n1:) ; n2 = index(buffer,'(') - 1
      read(buffer(1:n2),*) alpha
      n1 = index(buffer,'=') + 1 ; buffer = buffer(n1:) ; n2 = index(buffer,'(') - 1
      read(buffer(1:n2),*) beta
      n1 = index(buffer,'=') + 1 ; buffer = buffer(n1:) ; n2 = index(buffer,'(') - 1
      read(buffer(1:n2),*) gamma
    end if

! Now extract the fractional coordinates of the basis
    natoms = 0
    buffer = ''
    do while (index(buffer,'Name')==0)
      read(ic,'(a)') buffer
    end do
    read(ic,'(a)') buffer
    do while ((len(trim(adjustl(buffer)))/=0).and.(buffer(1:1)/='1'))
      natoms = natoms + 1
      if (natoms==1) then
        allocate(caxyz(natoms),occ(natoms),xf(natoms),yf(natoms),zf(natoms),atom_name(natoms))
      else
        allocate(caxyzt(natoms-1),occt(natoms-1),xft(natoms-1),yft(natoms-1),zft(natoms-1),atom_namet(natoms-1))     ! Here we create a temporary array ...
        caxyzt = caxyz                 ! ... into which we copy the contents of cxyz.
        occt = occ
        atom_namet = atom_name
        xft = xf
        yft = yf
        zft = zf
        if (allocated(caxyz)) deallocate(caxyz,occ,xf,yf,zf,atom_name)   ! Now we deallocate cxyz ...
        allocate(caxyz(natoms),occ(natoms),xf(natoms),yf(natoms),zf(natoms),atom_name(natoms))        ! ... so that we can reallocate it one item longer.
        caxyz(1:(natoms-1)) = caxyzt   ! Finally we copy the contents back into cxyz ...
        occ(1:(natoms-1)) = occt
        atom_name(1:(natoms-1)) = atom_namet
        xf(1:(natoms-1)) = xft
        yf(1:(natoms-1)) = yft
        zf(1:(natoms-1)) = zft
        if (allocated(caxyzt)) deallocate(caxyzt,occt,xft,yft,zft,atom_namet) ! ... and deallocate the temporary array.
      end if
!      do i = 1,8
!        if ( (ichar(buffer(i:i))>=ichar('0')) .and. (ichar(buffer(i:i))<=ichar('9')) ) buffer(i:i) = ' '
!      end do

      ! Revised algorithm to deal with atom labels ending in characters
      buffer = adjustl(buffer)
      i = index(buffer, ' ')
      do j = 1, i-1
         if ( ichar(buffer(j:j)) >= ichar('0') .and. &
              ichar(buffer(j:j)) <= ichar('9') ) & ! we've got to a number!
              buffer(j:i) = ' ' ! erase the rest of the first "word"
      end do

!WS Correction for atom labels in capital leters i.e. BA -> Ba
      if (ichar(buffer(2:2)) >= ichar('A').and.ichar(buffer(2:2)) <= ichar('Z') ) then
          buffer(2:2) = achar(ichar(buffer(2:2)) + 32)
      end if

!     Reading information about occupancy      
      if (allocated(citem)) deallocate(citem)
      call remove_errors(buffer)
      call remove_spaces_and_other(buffer)

!Yuanpeng - The codes below waw revised to account for the situation where the 
!site symbol in TBL file is given as something like "-42M 001".	  
	  if (allocated(citem_temp)) deallocate(citem_temp)
	  allocate(citem_temp(11))
	  read(buffer,*,iostat=ios) citem_temp
	  
	  if (citem_temp(7)(1:1)=="0") then
		citem_temp(6) = citem_temp(6)//citem_temp(7)
		citem_temp(7) = citem_temp(8)
		citem_temp(8) = citem_temp(9)
		citem_temp(9) = citem_temp(10)
		citem_temp(10) = citem_temp(11)
	  endif

	  allocate(citem(10))
	  
	  do i = 1,10
		citem(i) = citem_temp(i)
	  enddo
	  
      !read(buffer,*,iostat=ios) citem
! Yuanpeng's editing finished here.
	  read(citem(1),*) atom_name(natoms)
      read(citem(2),*) xf(natoms)
      read(citem(3),*) yf(natoms)
      read(citem(4),*) zf(natoms)
      read(citem(10),*) occ(natoms)
      !occ(natoms) = real(trim(citem(10)))
      
      caxyz(natoms) = adjustl(buffer(1:max_atom_string))
      read(ic,'(a)') buffer
    end do

    call split_occupancy

    do i = 1,natoms
      call remove_errors(caxyz(i))
      call remove_spaces(caxyz(i))
      ! >>>>>>>>>>>>>> Yuanpeng here -> the following line
      ! was removed since when reading in the TBL file (e.g the one
      ! given in tutorial 4), this line will cause problems. >>>>>> 
      !call remove_stuff(caxyz(i))
      ! <<<<<<<<<<<<<< Yuanpeng's editing finishes here. <<<<<<<<<<
      caxyz(i) = caxyz(i)(1:2)//'     '//caxyz(i)(3:max_atom_string-5)
    end do

    if (ldiag) then
      write(main_output,'(a)') 'Output from get_structure_tbl'
      write(main_output,'(a/)') '============================='
      write(main_output,'(a)') 'Crystal system is '//trim(system)
      write(main_output,'(a,3f8.4)') 'a,b,c = ',a,b,c
      write(main_output,'(a,3f9.4)') 'alpha,beta,gamma = ',alpha,beta,gamma
      write(main_output,'(a,i0)') 'Number of atoms = ',natoms
      do i = 1,natoms
        write(main_output,'(a)') trim(caxyz(i))
      end do
      write(main_output,'(/)')
    end if

    if (allocated(citem)) deallocate(citem)
    if (allocated(occ)) deallocate(occ,xf,yf,zf,atom_name)

    return

  end subroutine get_structure_tbl

!===============================================================================

  subroutine get_exp_spins
! ==============================

!--------------------------------------------------------------------------
!
!     Parses the structure information from a GSAS EXP file
!     paying particular attention to the magnetic spins
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use channel_numbers

    implicit none

    integer :: i,n,ierror,iatom
    character(len=1) :: ctype
    character(len=132) :: buffer

    open(imag,file=trim(filemag),status='old')
    main_loop: do
    read(imag,'(a)',iostat=ierror) buffer
    if (ierror>0) then
      write(6,'(a)') 'Problem locating lattice parameters from the exp file'
      stop
    end if
    if ((index(buffer,'CRS1')/=0).and.(index(buffer,'NATOM')/=0)) then
      n = index(buffer,' ',.true.)
      read(buffer(n:len_trim(buffer)),*) natoms
      allocate(caxyz(natoms),cspin(natoms))
      caxyz = '' ; cspin = ''
    else if ((index(buffer,'CRS1')/=0).and.(index(buffer,'ABC')/=0)) then
      n = index(buffer,'ABC') + len('ABC')
      read(buffer(n:),*) a,b,c
    else if ((index(buffer,'CRS1')/=0).and.(index(buffer,'ANGLES')/=0)) then
      n = index(buffer,'ANGLES') + len('ANGLES')
      read(buffer(n:),*) alpha,beta,gamma
      exit main_loop
    else if ((index(buffer,'CRS1')/=0).and.(index(buffer,'AT')/=0)) then
      n = index(buffer,'AT') + len('AT')
      read(buffer(n:n+3),*) iatom
      ctype = buffer(n+4:n+4)
      if (ctype=='A') caxyz(iatom) = buffer(15:52)
      if (ctype=='M') cspin(iatom) = buffer(15:44)
    end if
    end do main_loop

    if (ldiag) then
    write(main_output,'(/a/)') '=== Output from get_exp_spins ==='
    do i = 1,natoms
      if (len(cspin(i))>0) then
        write(main_output,'(a)') caxyz(iatom)//' M: '//cspin(i)
      else
        write(main_output,'(a)') caxyz(iatom)
      end if
    end do
    end if

    close(imag)

    return

  end subroutine get_exp_spins

!===============================================================================

  subroutine get_structure_cif
! ============================

!--------------------------------------------------------------------------
!
!     Parses the structure information from a CIF file
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use annotations
    use cif_data
    use cif_stuff
    use channel_numbers

    implicit none

! >>>>>>>>>>>>>>>>>>>>> Yuanpeng here -> Interface for the function >>>>>>>>>>>>
! determining whether the given string is numeric or not. >>>>>>>>>>>>>>>>>>>>>>
    INTERFACE 
        FUNCTION is_numeric(string)
            CHARACTER(len=*), INTENT(IN) :: string
            LOGICAL :: is_numeric
            REAL :: x
            INTEGER :: e
        END FUNCTION
    END INTERFACE
! <<<<<<<<<<<<<<<<<<<<< Yuanpeng finishes here <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    integer :: i,j,ios,nlines,nloops,n,natomloop,natomst,nsymloop,nld
    integer :: i_atom_site_label,i_atom_site_type_symbol,i_atom_site_fract_x, &
             i_atom_site_fract_y,i_atom_site_fract_z,i_atom_site_Cartn_x, &
             i_atom_site_Cartn_y,i_atom_site_Cartn_z,i_atom_site_symmetry_multiplicity, &
             i_atom_site_Wyckoff_symbol,i_atom_site_B_iso_or_equiv, &
             i_atom_site_occupancy
    integer :: i_symmetry_equiv_pos_site_id,i_space_group_symop_id, &
             i_symmetry_equiv_pos_as_xyz,i_space_group_symop_operation_xyz

    integer, allocatable :: loopstart(:),looplines(:),loopdescriptors(:)
    integer :: n_end_line,n_start_line
    integer, parameter :: buffer_length = 80
    character(len=buffer_length), allocatable :: line(:)
    character(len=buffer_length) :: buffer
    character(len=132), allocatable :: catoms(:)
    character(len=24), allocatable :: citem(:)
    logical :: lindescriptors

    i_atom_site_label = 0
    i_atom_site_type_symbol = 0
    i_atom_site_fract_x = 0
    i_atom_site_fract_y = 0
    i_atom_site_fract_z = 0
    i_atom_site_Cartn_x = 0
    i_atom_site_Cartn_y = 0
    i_atom_site_Cartn_z = 0
    i_atom_site_symmetry_multiplicity = 0
    i_atom_site_Wyckoff_symbol = 0
    i_atom_site_B_iso_or_equiv = 0
    i_atom_site_occupancy = 0
    i_symmetry_equiv_pos_site_id = 0
    i_space_group_symop_id = 0
    i_symmetry_equiv_pos_as_xyz = 0
    i_space_group_symop_operation_xyz = 0

! First count number of lines
    nlines = 0
    ios = 0
    line_count: do while (ios==0)
    read(ic,'(a)',iostat=ios) buffer
    if (ios/=0) exit line_count
    if (len(trim(buffer))==0) cycle line_count
    buffer = adjustl(buffer)
    if (buffer(1:1)=='#') cycle line_count
    nlines = nlines + 1
    if ((buffer(1:5)=='loop_').and.(len(trim(buffer))>5)) then
      buffer = adjustl(buffer(6:))
      do while (buffer(1:1)=='_')
        nlines = nlines + 1
        j = index(buffer,' ')
        buffer = adjustl(buffer(j:))
      end do
    end if
    end do line_count

! Put all the lines of the CIF file into the lines array
    allocate(line(nlines))
    rewind(ic)
    nlines = 0
    ios = 0
    line_read: do while (ios==0)
    read(ic,'(a)',iostat=ios) buffer
    if (ios/=0) exit line_read
    buffer = adjustl(buffer)
    if (len(trim(buffer))==0) cycle line_read
    if (buffer(1:1)=='#') cycle line_read
    nlines = nlines + 1
    line(nlines) = adjustl(buffer)
    if ((buffer(1:5)=='loop_').and.(len(trim(buffer))>5)) then
      line(nlines) = buffer(1:5)
      buffer = adjustl(buffer(6:))
      do while (buffer(1:1)=='_')
        nlines = nlines + 1
        j = index(buffer,' ')
        line(nlines) = adjustl(buffer(1:j-1))
        buffer = adjustl(buffer(j:))
      end do
    end if
    end do line_read

! Remove any tab characters
    do i = 1,nlines
    buffer = line(i)
    do while (index(buffer,char(9))>0)
      j = index(buffer,char(9))
      buffer(j:j) = ' '
    end do
    line(i) = adjustl(buffer)
    end do

! Count the number of loop_ lines in the CIF file
    nloops = 0
    do i = 1,nlines
    if (trim(line(i))=='loop_') nloops = nloops + 1
    end do
    allocate(loopstart(nloops),looplines(nloops),loopdescriptors(nloops))

! Take note of the number of each line that is a loop_ line
    nloops = 0
    do i = 1,nlines
    if (trim(line(i))=='loop_') then
      nloops = nloops + 1
      loopstart(nloops) = i
    end if
    end do

! Count number of lines within each loop_, and number of descriptor lines
    looplines = 0 ; loopdescriptors = 0
    lindescriptors = .true.
    do i = 1,nloops
    lindescriptors = .true.
    lines_per_loop: do j = loopstart(i)+1,nlines
      buffer = trim(line(j))
      if ((buffer(1:1)=='_').and.lindescriptors) then
        loopdescriptors(i) = loopdescriptors(i) + 1
        looplines(i) = looplines(i) + 1
        cycle lines_per_loop
      else
        lindescriptors = .false.
      end if
      if ((index(buffer,'loop_')>0).and..not.lindescriptors) exit lines_per_loop
      looplines(i) = looplines(i) + 1
    end do lines_per_loop
    end do

! Look through line list for lattice parameters and space group
    do i = 1,nlines
    buffer = line(i)
    if (buffer(1:14)=='_cell_length_a') then
      buffer = adjustl(buffer(15:))
      n = index(buffer,'(') - 1
      if (n>0) then
        read(buffer(1:n),*) a
      else
        read(buffer,*) a
      end if
    end if
    if (buffer(1:14)=='_cell_length_b') then
      buffer = adjustl(buffer(15:))
      n = index(buffer,'(') - 1
      if (n>0) then
        read(buffer(1:n),*) b
      else
        read(buffer,*) b
      end if
    end if
    if (buffer(1:14)=='_cell_length_c') then
      buffer = adjustl(buffer(15:))
      n = index(buffer,'(') - 1
      if (n>0) then
        read(buffer(1:n),*) c
      else
        read(buffer,*) c
      end if
    end if
    if (buffer(1:17)=='_cell_angle_alpha') then
      buffer = adjustl(buffer(18:))
      n = index(buffer,'(') - 1
      if (n>0) then
        read(buffer(1:n),*) alpha
      else
        read(buffer,*) alpha
      end if
    end if
    if (buffer(1:16)=='_cell_angle_beta') then
      buffer = adjustl(buffer(17:))
      n = index(buffer,'(') - 1
      if (n>0) then
        read(buffer(1:n),*) beta
      else
        read(buffer,*) beta
      end if
    end if
    if (buffer(1:17)=='_cell_angle_gamma') then
      buffer = adjustl(buffer(18:))
      n = index(buffer,'(') - 1
      if (n>0) then
        read(buffer(1:n),*) gamma
      else
        read(buffer,*) gamma
      end if
    end if
    if (buffer(1:30)=='_symmetry_space_group_name_H-M') then
      buffer = adjustl(buffer(31:))
      metadata_spaceg = trim(buffer)
      if (metadata_spaceg(1:1)=='''') metadata_spaceg = metadata_spaceg(2:)
      n = index(metadata_spaceg,'''')
      if (n>0) metadata_spaceg = metadata_spaceg(1:n-1)
    end if
    if (buffer(1:25)=='_space_group_name_H-M_alt') then
      buffer = adjustl(buffer(26:))
      metadata_spaceg = trim(buffer)
      if (metadata_spaceg(1:1)=='''') metadata_spaceg = metadata_spaceg(2:)
      n = index(metadata_spaceg,'''')
      if (n>0) metadata_spaceg = metadata_spaceg(1:n-1)
    end if
    end do

! Scan through lines for some of the metadata
    do i = 1,nlines
    buffer = line(i)
    if (trim(buffer)=='_chemical_name_systematic') then
      metadata_title = adjustl(trim(line(i+1)))
      if (metadata_title(1:1)=='''') metadata_title = metadata_title(2:)
      n = index(metadata_title,'''')
      if (n>0) metadata_title = metadata_title(1:n-1)
    end if
    if (trim(buffer)=='_chemical_formula_structural') then
      metadata_material = trim(line(i+1))
      if (metadata_material(1:1)=='''') metadata_material = metadata_material(2:)
      n = index(metadata_material,'''')
      if (n>0) metadata_material = metadata_material(1:n-1)
    end if
    end do

! Look for information about author
   metadata_author = ''
   authorloop: do i = 1,nloops
     buffer = line(loopstart(i)+1)
     if (buffer(1:12)/='_publ_author') cycle authorloop
     if ((looplines(i)-loopdescriptors(i))==1) then
     metadata_author = trim(line(loopstart(i)+2))
     exit authorloop
     end if
     if ((looplines(i)-loopdescriptors(i))==2) then
     if (trim(buffer)=='_publ_author_name') then
       metadata_author = trim(line(loopstart(i)+3))
       exit authorloop
     else
       metadata_author = trim(line(loopstart(i)+looplines(i)))
       exit authorloop
     end if
     end if
   end do authorloop
   if (ldiag) then
     write(main_output,*)
     write(main_output,'(a)') 'Output from get_structure_cif'
     do i = 1,nlines
     write(main_output,'(a)') trim(line(i))
     end do
     write(main_output,*)
     write(main_output,'(a,i0)') 'nloops: ',nloops
     do i = 1,nloops
     write(main_output,*) i,loopstart(i),loopdescriptors(i),looplines(i)
     end do
   end if

! Look for information about symmetry
! Programming note. This is known to work for the case where the CIF file has one
! operator per line, and with a number in front of the operator, eg lines of the form
! 1 'x,y,z'
! It is possible to put more than one operator per line, and this is coded in, but not
! really tested. Moreover, any other cases may not have been anticipated properly.
   nsym = 0 ; nsymloop = 0
   symmetryloop: do i = 1,nloops
     do j = loopstart(i)+1,loopstart(i)+loopdescriptors(i)
     if (index(line(j),'_symmetry_equiv_pos_')>0) nsymloop = i
     if (index(line(j),'_space_group_symop_')>0) nsymloop = i
     end do
   end do symmetryloop
   if (nsymloop>0) then
     nsym = looplines(nsymloop) - loopdescriptors(nsymloop)
     allocate(coperators(nsym))
     i_symmetry_equiv_pos_site_id = 0
     i_space_group_symop_id = 0
     i_symmetry_equiv_pos_as_xyz = 0
     i_space_group_symop_operation_xyz = 0
     do i = loopstart(nsymloop)+1,loopstart(nsymloop)+loopdescriptors(nsymloop)
       buffer = adjustl(line(i))
       if (trim(buffer)=='_symmetry_equiv_pos_site_id') &
             i_symmetry_equiv_pos_site_id = i - loopstart(nsymloop)
       if (trim(buffer)=='_space_group_symop_id') &
             i_space_group_symop_id = i - loopstart(nsymloop)
       if (trim(buffer)=='_symmetry_equiv_pos_as_xyz') &
              i_symmetry_equiv_pos_as_xyz = i - loopstart(nsymloop)
       if (trim(buffer)=='_space_group_symop_operation_xyz') &
              i_space_group_symop_operation_xyz = i - loopstart(nsymloop)
     end do
     nld = loopdescriptors(nsymloop)
     if (allocated(citem)) deallocate(citem)
     allocate(citem(nld))
     do j = 1,nsym
       buffer = trim(line(loopstart(nsymloop)+loopdescriptors(nsymloop)+j))
       buffer = adjustl(buffer)
       call remove_space(buffer)
       ! call remove_space(buffer) 

       ! The above line was added back by Yuanpeng&Martin since we
       ! we found without this line, reading in some CIF file 
       ! generated by VESTA may encounter some problems.
       ! Not sure whether bringing the line back is a good idea or not.
       
       ! The above line removed by AEP - it doesn't play nicely with
       ! the following loop that actually looks for spaces. I wonder
       ! if it was introduced in Nov 2009?!
       
       do i = 1,nld   !!! I am not sure what this loop is doing, Nov 2009
!         n = index(buffer,' ')
         citem(i) = trim(buffer)
!         buffer = adjustl(buffer(n:))
       end do
       if (i_symmetry_equiv_pos_as_xyz>0) buffer = adjustl(citem(i_symmetry_equiv_pos_as_xyz))
       if (i_space_group_symop_id>0)  buffer = adjustl(citem(i_space_group_symop_id))
       if (i_symmetry_equiv_pos_as_xyz>0)  buffer = adjustl(citem(i_symmetry_equiv_pos_as_xyz))
       if (i_space_group_symop_operation_xyz>0)  buffer = adjustl(citem(i_space_group_symop_operation_xyz))
       !if (buffer(1:1)=='''') buffer = buffer(2:)
       ! >>>>>>>>>>>>>>> The line above was changed to the following
       ! two lines by Yuanpeng. For some CIF files, the symmetry operator
       ! line starts with an index number, which causes problem when reading
       ! in the operators >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       n = index(buffer,'''')
       if (n>0) buffer = buffer(n+1:)
       n = index(buffer,'''')
       if (n>0) buffer = buffer(1:n-1)
        ! >>>>>>>>>> Yuanpeng here -> Sometimes the symmetry operator is given
        ! in this way -> 16 Z,X,-Y without single quote sign. The following
        ! code is for dealing with such situation. >>>>>>>>>>>>>>>>>>>>>>>>>>>
       if (.NOT. buffer(2:2)/="//") then
          write(*,*) buffer(2:2)
          i=0
          do while (is_numeric(buffer(i+1:i+1)))
           i = i + 1
          end do
          buffer = buffer(i+1:)
       end if
       ! <<<<<<<<<<< Yuanpeng finishes here <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
       coperators(j) = trim(buffer)
     end do
   end if
!     end if
!     if (loopdescriptors(i)==1) then   ! Here we assume lines lines containing any number
                                     ! of operators
!       nsymlines = looplines(i) - loopdescriptors(i)
!       nsym = 0
!       do j = 1,nsymlines
!         buffer = adjustl(trim(line(loopstart(i)+loopdescriptors(i)+j)))
!         do while (len(trim(buffer))>0)
!           n = index(buffer,' ')
!           nsym = nsym + 1
!           csymoplines(nsym) = buffer(1:n-1)
!           buffer = adjustl(buffer(n:))
!         end do
!       end do
!       allocate(coperators(nsym))
!       coperators = csymoplines(1:nsym)
!       do j = 1,nsym
!         buffer = adjustl(coperators(j))
!         if (buffer(1:1)=='''') buffer = buffer(2:)
!         n = index(buffer,'''')
!         if (n>0) buffer = buffer(1:n-1)
!         coperators(j) = buffer
!       end do
!     end if
!   end do

! Finally get the atom list
   i_atom_site_label = 0
   i_atom_site_type_symbol = 0
   i_atom_site_fract_x = 0
   i_atom_site_fract_y = 0
   i_atom_site_fract_z = 0
   i_atom_site_Cartn_x = 0
   i_atom_site_Cartn_y = 0
   i_atom_site_Cartn_z = 0
   atomloop: do i = 1,nloops
     do j = loopstart(i)+1,loopstart(i)+loopdescriptors(i)
       if (trim(line(j))=='_atom_site_fract_x') then
         natomloop = i
         exit atomloop
       end if
     end do
   end do atomloop
!    write(6,'(a,i0)') 'natomloop = ',natomloop
!    write(6,'(a,i0)') 'loopstart(natomloop) = ',loopstart(natomloop)
!    write(6,'(a,i0)') 'loopdescriptors(natomloop) = ',loopdescriptors(natomloop)
!    write(6,'(a,i0)') 'looplines(natomloop) = ',looplines(natomloop)
   
   do i = loopstart(natomloop)+1,loopstart(natomloop)+loopdescriptors(natomloop)
     buffer = line(i)
     if (trim(buffer)=='_atom_site_label') i_atom_site_label = i - loopstart(natomloop)
     if (trim(buffer)=='_atom_site_type_symbol') &
           i_atom_site_type_symbol = i - loopstart(natomloop)
     if (trim(buffer)=='_atom_site_symmetry_multiplicity') &
           i_atom_site_symmetry_multiplicity = i - loopstart(natomloop)
     if (trim(buffer)=='_atom_site_Wyckoff_symbol') &
           i_atom_site_Wyckoff_symbol = i - loopstart(natomloop)
     if (trim(buffer)=='_atom_site_fract_x') i_atom_site_fract_x = i - loopstart(natomloop)
     if (trim(buffer)=='_atom_site_fract_y') i_atom_site_fract_y = i - loopstart(natomloop)
     if (trim(buffer)=='_atom_site_fract_z') i_atom_site_fract_z = i - loopstart(natomloop)
     if (trim(buffer)=='_atom_site_Cartn_x') i_atom_site_Cartn_x = i - loopstart(natomloop)
     if (trim(buffer)=='_atom_site_Cartn_y') i_atom_site_Cartn_y = i - loopstart(natomloop)
     if (trim(buffer)=='_atom_site_Cartn_z') i_atom_site_Cartn_z = i - loopstart(natomloop)
     if (trim(buffer)=='_atom_site_B_iso_or_equiv') &
           i_atom_site_B_iso_or_equiv = i - loopstart(natomloop)
     if (trim(buffer)=='_atom_site_occupancy') &
           i_atom_site_occupancy = i - loopstart(natomloop)
   end do
   n = loopdescriptors(natomloop)
   nlines = looplines(natomloop) - loopdescriptors(natomloop)
   n_start_line = loopstart(natomloop)+loopdescriptors(natomloop)+1
   n_end_line = n_start_line + nlines - 1
! write(6,'(a,i0)') 'nlines = ',nlines
! write(6,'(a,i0)') 'n_start_line = ',n_start_line
! write(6,'(a,i0)') 'n_end_line   = ',n_end_line
   call count_cif_items(line(n_start_line:n_end_line),n,natomst)
   if (allocated(citem)) deallocate(citem)
   allocate(citem(n))
   natoms = natomst
!!   natoms = looplines(natomloop) - loopdescriptors(natomloop)
   allocate(catoms(natoms))
   allocate(xf(natoms),yf(natoms),zf(natoms))
   allocate(xo(natoms),yo(natoms),zo(natoms))
   allocate(occ(natoms))
   allocate(atom_name(natoms))
!    write(6,'(a,i0)') 'Number of atoms = ',natoms
   do i = 1,natoms
!!     j = loopstart(natomloop)+loopdescriptors(natomloop)+i
     buffer = atom_data_lines(i)(1:buffer_length)
     read(buffer,*,iostat=ios) citem
!     write(6,*) ios
!write(6,'(a)') trim(adjustl(citem(i_atom_site_label)))
!write(6,'(a)') trim(adjustl(citem(i_atom_site_type_symbol)))
!stop
     if (i_atom_site_label>0) atom_name(i) = trim(adjustl(citem(i_atom_site_label)))
     if (i_atom_site_type_symbol>0) atom_name(i) = trim(adjustl(citem(i_atom_site_type_symbol)))
     if (i_atom_site_fract_x>0) then
     buffer = adjustl(trim(citem(i_atom_site_fract_x)))
     n = index(buffer,'(')
     if (n>0) buffer = buffer(1:n-1)
     read(buffer,*) xf(i)
     end if
     if (i_atom_site_fract_y>0) then
     buffer = adjustl(trim(citem(i_atom_site_fract_y)))
     n = index(buffer,'(')
     if (n>0) buffer = buffer(1:n-1)
     read(buffer,*) yf(i)
     end if
     if (i_atom_site_fract_z>0) then
     buffer = adjustl(trim(citem(i_atom_site_fract_z)))
     n = index(buffer,'(')
     if (n>0) buffer = buffer(1:n-1)
     read(buffer,*) zf(i)
     end if
     if (i_atom_site_Cartn_x>0) then
     buffer = adjustl(trim(citem(i_atom_site_Cartn_x)))
     n = index(buffer,'(')
     if (n>0) buffer = buffer(1:n-1)
     read(buffer,*) xo(i)
     end if
     if (i_atom_site_Cartn_y>0) then
     buffer = adjustl(trim(citem(i_atom_site_Cartn_y)))
     n = index(buffer,'(')
     if (n>0) buffer = buffer(1:n-1)
     read(buffer,*) yo(i)
     end if
     if (i_atom_site_Cartn_z>0) then
     buffer = adjustl(trim(citem(i_atom_site_Cartn_z)))
     n = index(buffer,'(')
     if (n>0) buffer = buffer(1:n-1)
     read(buffer,*) zo(i)
     end if
     if (i_atom_site_occupancy>0) then
        buffer = adjustl(trim(citem(i_atom_site_occupancy)))
        n = index(buffer,'(')
        if (n>0) buffer = buffer(1:n-1)
        read(buffer,*) occ(i)
     else
        occ(i) = 1.0 
     end if
   end do
   
! WS extra code

  if (any(occ(:) < 1.0)) then
      call split_occupancy
  endif
  

! endo of WS
   
   
   
   allocate(caxyz(natoms))
   do i = 1,natoms
     write(caxyz(i),'(a8,3f12.5)') adjustl(trim(atom_name(i))),xf(i),yf(i),zf(i)
   end do

   if (ldiag) then
     write(main_output,'(a)') 'Output from get_structure_cif'
     write(main_output,'(a,i0)') 'Number of atoms = ',natoms
     write(main_output,'(a)') 'Atomic coordinates'
     do i = 1,natoms
     write(main_output,'(/a)') trim(caxyz(i))
     end do
     write(main_output,'(a)') 'Symmetry operators'
     do i = 1,nsym
     write(main_output,'(/a)') coperators(i)
     end do
   end if
!   write(6,'(a,i0)') 'Number of space group operators = ',nsym
   lattice_info(1) = nsym
   lattice_info(2) = 0     ! Although there may be a centre of symmetry, this will
                         ! be in the list of symmetry operators
   centre = 'P'            ! May not be, but again, centre is contained with the
                         ! list of symmetry operators
   if (allocated(xf)) deallocate(xf,yf,zf)
   if (allocated(xo)) deallocate(xo,yo,zo)
   if (allocated(occ)) deallocate(occ)
   if (allocated(atom_name)) deallocate(atom_name)
   if (allocated(atom_label)) deallocate(atom_label)
   if (allocated(catom_label)) deallocate(catom_label)

  end subroutine get_structure_cif

!===============================================================================

  subroutine get_structure_sff
! ============================

!--------------------------------------------------------------------------
!
!     Parses the structure information from a simple file format (sff) file
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use channel_numbers

    implicit none

    integer :: i,n
    character(len=80) :: buffer
    character(len=8) :: catom

! Read cell parameters
    buffer = ''
    do while (index(buffer,'cell')==0)
      read(ic,'(a)') buffer
    end do
    read(ic,*) a,b,c,alpha,beta,gamma
    rewind(ic)

    buffer = ''
    do while (index(buffer,'symmetry')==0)
      read(ic,'(a)') buffer
    end do
    buffer= adjustl(buffer)
    read(buffer(9:),*) nsym
    allocate(coperators(nsym))
    lattice_info(1) = nsym ; lattice_info(2) = 0
    do i = 1,nsym
      read(ic,'(a)') coperators(i)
    end do
    rewind(ic)

    buffer = ''
    do while (index(buffer,'atoms')==0)
      read(ic,'(a)') buffer
    end do
    buffer= adjustl(buffer)
    read(buffer(6:),*) natoms
    allocate(caxyz(natoms))
    do i = 1,natoms
      read(ic,'(a)') buffer
      buffer = adjustl(buffer)
      n = index(buffer,' ')
      catom = buffer(1:n-1)
      buffer = catom//buffer(n:)
      caxyz(i) = buffer(1:max_atom_string)
    end do

   if (ldiag) then
     write(main_output,*)
     write(main_output,'(a)') 'Output from get_structure_sff'
     write(main_output,'(a,i0)') 'Number of atoms = ',natoms
     do i = 1,natoms
     write(main_output,'(a)') trim(caxyz(i))
     end do
     do i = 1,nsym
     write(main_output,'(a)') coperators(i)
     end do
   end if

  end subroutine get_structure_sff

!===============================================================================

  subroutine get_structure_resfile
! ================================

!--------------------------------------------------------------------------
!
!     Parses the structure information from a SHELX res format file
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use channel_numbers

    implicit none

    integer :: i,n,iatom,ierror
    character(len=132) :: buffer,csfac
    character(len=8) :: catom
    double precision :: temp(4),wavelength

    nsym = 1

    readfile_general: do
      read(ic,'(a)',iostat=ierror) buffer
      if (ierror/=0) exit readfile_general
      buffer = adjustl(buffer)
      n = index(buffer,' ') - 1
      select case(buffer(1:n))
        case('CELL')
          read(buffer(n+1:),*) wavelength,a,b,c,alpha,beta,gamma
        case('LATT')
          read(buffer(n+1:),*) i
          if (i<0) then
            lattice_info(2) = 0
          else
            lattice_info(2) = 1
          end if
        case('SYMM')
          nsym = nsym + 1
        case('SFAC')
          csfac = adjustl(buffer(n+1:))
        case default
          continue
      end select
    end do readfile_general

    rewind(ic)
    
    allocate(coperators(nsym))

! Capture the symmetry information
    i = 1
    lattice_info(1) = nsym
    coperators(1) = '+X,+Y,+Z'
    readfile_symmetry: do
      read(ic,'(a)',iostat=ierror) buffer
      if (ierror/=0) exit readfile_symmetry
      buffer = adjustl(buffer)
      n = index(buffer,' ') - 1
      select case(buffer(1:n))
        case('SYMM')
          i = i + 1
          read(buffer(n+1:),'(a)') coperators(i)
        case default
          continue
      end select
    end do readfile_symmetry
    rewind(ic)
    
! Now find the number of atoms
   natoms = 0
   ierror = 0
   readfile_atoms1: do
      read(ic,'(a)',iostat=ierror) buffer
      if (ierror/=0) exit readfile_atoms1
      buffer = adjustl(buffer)
      if(index(buffer,',')>0) cycle readfile_atoms1
      if(index(buffer,'/')>0) cycle readfile_atoms1
      n = index(buffer,' ')
      ierror = 0
      read(buffer(n:),*,iostat=ierror) i,temp
      if (ierror/=0) cycle readfile_atoms1
      natoms = natoms + 1
   end do readfile_atoms1

  allocate(caxyz(natoms))    

! Now find the actual atoms
   rewind(ic)
   iatom = 0
   readfile_atoms2: do
      read(ic,'(a)',iostat=ierror) buffer
      if (ierror/=0) exit readfile_atoms2
      buffer = adjustl(buffer)
      if(index(buffer,',')>0) cycle readfile_atoms2
      if(index(buffer,'/')>0) cycle readfile_atoms2
      n = index(buffer,' ')
      read(buffer(n:),*,iostat=ierror) i,temp
      if (ierror/=0) cycle readfile_atoms2
      iatom = iatom + 1
      catom = buffer(1:2)
      if(ichar(catom(2:2))>=ichar('1').and.ichar(catom(2:2))<=ichar('9')) catom(2:2) = ' '
      write(buffer,'(3f12.5)') temp(1:3)
      caxyz(iatom) = catom//trim(buffer)
   end do readfile_atoms2 


   if (ldiag) then
     write(main_output,*)
     write(main_output,'(a)') 'Output from get_structure_resfile'
     write(main_output,'(a,i0)') 'Number of atoms = ',natoms
     do i = 1,natoms
       write(main_output,'(a)') trim(caxyz(i))
     end do
     do i = 1,nsym
       write(main_output,'(a)') coperators(i)
     end do
   end if

   return

   end subroutine get_structure_resfile

!===============================================================================


  subroutine label_atoms
! ======================

!--------------------------------------------------------------------------
!
!     This subroutine generates the atom labels, eg Na1 and Na2, and then
!     O1, O2 and O3
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use utilities
    use channel_numbers

    implicit none
    
    integer,parameter :: ntext = 8
    integer :: i,j,k,m,n,ierror
    character(len=ntext) :: catominfo
    character(len=(ntext-2)) :: clabel
    character(len=2) :: ci,cj,cel
    integer, allocatable :: nlabel(:)
    
    allocate(nlabel(natoms))

    do i = 1,natoms
      catominfo = adjustl(caxyz(i)(1:ntext))
      call extract_atom_name_label(catominfo,cel,clabel)
      caxyz(i)(1:ntext) = cel//clabel
    end do

    if (lrelabel) then
      nlabel = 0
      do i = 1,natoms
        ci = caxyz(i)(1:2)
        n = 1
        do j = 1,i-1
          cj = caxyz(j)(1:2)
          if (ci==cj) n = n + 1
        end do
        nlabel(i) = n
      end do
      do i = 1,natoms
        write(caxyz(i)(3:ntext),'(i0)') nlabel(i)
      end do
    end if
    
!WS Apply labels based on the element type only    
    if (lelabel) then
      nlabel = 1
      do i = 1,natoms
        write(caxyz(i)(3:ntext),'(i0)') nlabel(i)
      end do
    end if
    
    if (luselabels) then
      do i = 1,natoms
        if (len_trim(caxyz(i)(3:ntext))==0) caxyz(i)(3:3) = '0'
        read(caxyz(i)(3:ntext),*,iostat=ierror) n
        if (ierror/=0) then
          write(6,'(a,i0)') 'Error detected reading numeric atom label for atom number ',i
          stop
        end if
      end do
    end if

   if (ldiag) then
     write(main_output,'(/a)') 'Output from label_atoms' 
     write(main_output,'(a/)') '=======================' 
     write(main_output,'(a,i0)') 'Number of atoms = ',natoms
     do i = 1,natoms
       write(main_output,'(i0,2x,a)') i,trim(caxyz(i))
     end do
     write(main_output,'(/)')
   end if

    deallocate(nlabel)

    return
    
  end subroutine label_atoms    

!===============================================================================

  subroutine generate_structure
! =============================

!--------------------------------------------------------------------------
!
!     This subroutine generates the crystal structure from the space
!     group operations and the input atomic coordinates.
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use utilities
    use channel_numbers

    implicit none

    double precision :: al,be,ga
    double precision :: s,volume,xr,yr,zr,x,y,z,xff,yff,zff,dxf,dyf,dzf, &
                        dxo,dyo,dzo,dr
    integer :: nmax,ix1,ix2,iy1,iy2,iz1,iz2,lattice,i,j,k,ij, &
               i2,niatoms,iel,ierror
    integer :: ifrac1(3),ifrac2(3)
    character(len=80) :: cbuffer
    character(len=8) :: can
    character(len=2) :: cx1,cx2,cy1,cy2,cz1,cz2
    character(len=2) :: c1(3),c2(3)
    character(len=26) :: cstring
    logical :: okay

    cell = 0.0
    n_elements = 0

! ... first form the lattice vectors
!     Comment: note the definition here
!     cell(i,j) has j define the vector, eg cell(:,1) points along the a lattice vector,
!     and i defines x,y,z, so that cell(3:2) is the z-component of the b lattice vector.
    al = alpha/57.295779512d0
    be = beta/57.295779512d0
    ga = gamma/57.295779512d0
    s = 0.5*(al+be+ga)
!    volume = 2.0*a*b*c*sqrt(sin(s)*sin(s-al)*sin(s-be)*sin(s-ga))
    volume = a*b*c*sqrt(1.0d0-cos(2.0d0*al)-cos(2.0d0*be)-cos(2.0d0*ga)) + 2.0d0*(cos(al)*cos(be)*cos(ga))
    cell(1,1) = volume/(sin(al)*b*c)
    cell(2,1) = a*cos(ga)*sin(al)
    cell(3,1) = a*cos(be)
    cell(1,2) = 0.0d0
    cell(2,2) = b*sin(al)
    cell(3,2) = b*cos(al)
    cell(1,3) = 0.0d0
    cell(2,3) = 0.0d0
    cell(3,3) = c

    if (.not.allocated(coperators)) call symmetry_error()  !Exit the program when no symmetry operators has been read

    nmax = 2*natoms*lattice_info(1)*(lattice_info(2)+1) ! Take account of centre of symmetry
                                                        ! doubling the number of atoms, and
                                                        ! multiply by 2 again in case we have
                                                        ! an A, B, C or I centre
    if (centre.eq.'P') nmax = nmax/2                    ! Divide by 2 if we have a P lattice
    if (centre.eq.'F') nmax = nmax*2                    ! and multiply in case of an F

    allocate(xf(nmax),yf(nmax),zf(nmax))
    allocate(xo(nmax),yo(nmax),zo(nmax))
    allocate(atom_name(nmax))
    allocate(atom_type(nmax))
!    if (luselabels) then
!      allocate(catom_label(nmax))
!    else
!      allocate(atom_label(nmax))
!    end if
    allocate(atom_label(nmax))
    allocate(reference_cell(nmax,3),reference_number(nmax))

    if (ldiag) then
      write(main_output,'(/a)') 'Output from generate_structure'
      write(main_output,'(a/)') '=============================='
      write(main_output,'(a,i0)') 'Number of atoms = ',natoms
!      write(main_output,'(a)') ' === Original configuration in fractional coordinates ==='
    end if
    
! ... read in the fractional coordinates, generate all the related
!     coordinates by looping over symmetry operations, the centre, and the
!     lattice types, and then check that there are no duplicates.
    niatoms = natoms ; natoms = 0
! write(6,'(a,i0,1x,i0)') 'Lattice: ',lattice_info
      do  i2 = 1,lattice_info(1)
        cstring = adjustl(coperators(i2))
! ... Parse the symmetry operator strings
        call convert_coperators(c1,c2,ifrac1,ifrac2,cstring)
        cx1 = c1(1) ; cy1 = c1(2); cz1 = c1(3)
        cx2 = c2(1) ; cy2 = c2(2); cz2 = c2(3)
        ix1 = ifrac1(1) ; iy1 = ifrac1(2) ; iz1 = ifrac1(3)
        ix2 = ifrac2(1) ; iy2 = ifrac2(2) ; iz2 = ifrac2(3)
        if (ldiag)  then
          write(main_output,'(a)') trim(cstring)
          write(main_output,'(5x,3(i1,a1,i1,2a2,2x))') &
          ix1,'/',ix2,cx1,cx2,iy1,'/',iy2,cy1,cy2,iz1,'/',iz2,cz1,cz2
        end if

!    if (ldiag) then
!      do i = 1,niatoms
!        cbuffer = adjustl(caxyz(i))
!        write(main_output,'(a)') trim(cbuffer)
!        write(main_output,'(a)') trim(cbuffer(13:))
!      end do
!    end if

    do i = 1,niatoms
      cbuffer = adjustl(caxyz(i))
!      n = index(cbuffer,' ')
      read(cbuffer(9:),*) xr,yr,zr
      j = nint(10000.0d0*xr)
      if (j==3333) xr = 1.0d0/3.0d0
      if ((j==6667).or.(j==6666)) xr = 2.0d0/3.0d0
      j = nint(10000.0d0*yr)
      if (j==3333) yr = 1.0d0/3.0d0
      if ((j==6667).or.(j==6666)) yr = 2.0d0/3.0d0
      j = nint(10000.0d0*zr)
      if (j==3333) zr = 1.0d0/3.0d0
      if ((j==6667).or.(j==6666)) zr = 2.0d0/3.0d0

      can = cbuffer(1:8)
!      if (ldiag) write(main_output,'(a,2x,3f10.4)') trim(can),xr,yr,zr

! ... Now apply the symmetry
        x = 0.0d0
        y = 0.0d0
        z = 0.0d0
        if (ix2.ne.0) x = dble(ix1)/dble(ix2)
        if (iy2.ne.0) y = dble(iy1)/dble(iy2)
        if (iz2.ne.0) z = dble(iz1)/dble(iz2)
        if ((cx1.eq.'+x').or.(cx1.eq.'+X')) x = x + xr
        if ((cx1.eq.' x').or.(cx1.eq.' X')) x = x + xr
        if ((cx1.eq.'-x').or.(cx1.eq.'-X')) x = x - xr
        if ((cx1.eq.'+y').or.(cx1.eq.'+Y')) x = x + yr
        if ((cx1.eq.' y').or.(cx1.eq.' Y')) x = x + yr
        if ((cx1.eq.'-y').or.(cx1.eq.'-Y')) x = x - yr
        if ((cx1.eq.'+z').or.(cx1.eq.'+Z')) x = x + zr
        if ((cx1.eq.' z').or.(cx1.eq.' Z')) x = x + zr
        if ((cx1.eq.'-z').or.(cx1.eq.'-Z')) x = x - zr
        if ((cy1.eq.'+x').or.(cy1.eq.'+X')) y = y + xr
        if ((cy1.eq.' x').or.(cy1.eq.' X')) y = y + xr
        if ((cy1.eq.'-x').or.(cy1.eq.'-X')) y = y - xr
        if ((cy1.eq.'+y').or.(cy1.eq.'+Y')) y = y + yr
        if ((cy1.eq.' y').or.(cy1.eq.' Y')) y = y + yr
        if ((cy1.eq.'-y').or.(cy1.eq.'-Y')) y = y - yr
        if ((cy1.eq.'+z').or.(cy1.eq.'+Z')) y = y + zr
        if ((cy1.eq.' z').or.(cy1.eq.' Z')) y = y + zr
        if ((cy1.eq.'-z').or.(cy1.eq.'-Z')) y = y - zr
        if ((cz1.eq.'+x').or.(cz1.eq.'+X')) z = z + xr
        if ((cz1.eq.' x').or.(cz1.eq.' X')) z = z + xr
        if ((cz1.eq.'-x').or.(cz1.eq.'-X')) z = z - xr
        if ((cz1.eq.'+y').or.(cz1.eq.' Y')) z = z + yr
        if ((cz1.eq.' y').or.(cz1.eq.'+Y')) z = z + yr
        if ((cz1.eq.'-y').or.(cz1.eq.'-Y')) z = z - yr
        if ((cz1.eq.'+z').or.(cz1.eq.'+Z')) z = z + zr
        if ((cz1.eq.' z').or.(cz1.eq.' Z')) z = z + zr
        if ((cz1.eq.'-z').or.(cz1.eq.'-Z')) z = z - zr
        if ((cx2.eq.'+x').or.(cx2.eq.'+X')) x = x + xr
        if ((cx2.eq.' x').or.(cx2.eq.' X')) x = x + xr
        if ((cx2.eq.'-x').or.(cx2.eq.'-X')) x = x - xr
        if ((cx2.eq.'+y').or.(cx2.eq.'+Y')) x = x + yr
        if ((cx2.eq.' y').or.(cx2.eq.' Y')) x = x + yr
        if ((cx2.eq.'-y').or.(cx2.eq.'-Y')) x = x - yr
        if ((cx2.eq.'+z').or.(cx2.eq.'+Z')) x = x + zr
        if ((cx2.eq.' z').or.(cx2.eq.' Z')) x = x + zr
        if ((cx2.eq.'-z').or.(cx2.eq.'-Z')) x = x - zr
        if ((cy2.eq.'+x').or.(cy2.eq.'+X')) y = y + xr
        if ((cy2.eq.' x').or.(cy2.eq.' X')) y = y + xr
        if ((cy2.eq.'-x').or.(cy2.eq.'-X')) y = y - xr
        if ((cy2.eq.'+y').or.(cy2.eq.'+Y')) y = y + yr
        if ((cy2.eq.' y').or.(cy2.eq.' Y')) y = y + yr
        if ((cy2.eq.'-y').or.(cy2.eq.'-Y')) y = y - yr
        if ((cy2.eq.'+z').or.(cy2.eq.'+Z')) y = y + zr
        if ((cy2.eq.' z').or.(cy2.eq.' Z')) y = y + zr
        if ((cy2.eq.'-z').or.(cy2.eq.'-Z')) y = y - zr
        if ((cz2.eq.'+x').or.(cz2.eq.'+X')) z = z + xr
        if ((cz2.eq.' x').or.(cz2.eq.' X')) z = z + xr
        if ((cz2.eq.'-x').or.(cz2.eq.'-X')) z = z - xr
        if ((cz2.eq.'+y').or.(cz2.eq.'+Y')) z = z + yr
        if ((cz2.eq.' y').or.(cz2.eq.' Y')) z = z + yr
        if ((cz2.eq.'-y').or.(cz2.eq.'-Y')) z = z - yr
        if ((cz2.eq.'+z').or.(cz2.eq.'+Z')) z = z + zr
        if ((cz2.eq.' z').or.(cz2.eq.' Z')) z = z + zr
        if ((cz2.eq.'-z').or.(cz2.eq.'-Z')) z = z - zr
        do j = 1,(lattice_info(2)+1)
          if (j.eq.2) then
            x = -x
            y = -y
            z = -z
          endif
          if (ladjust_xyz_signs) then
            do while (x<0.0d0)
              x = x + 1.0d0
            end do
            do while (y<0.0d0)
              y = y + 1.0d0
            end do
            do while (z<0.0d0)
              z = z + 1.0d0
            end do
            do while (x>=1.0d0)
              x = x - 1.0d0
            end do
            do while (y>=1.0d0)
              y = y - 1.0d0
            end do
            do while (z>=1.0d0)
              z = z - 1.0d0
            end do
          end if  
          lattice = 2   ! Correct for A, B, C and I
          if (centre.eq.'P') lattice = 1
          if (centre.eq.'F') lattice = 4
          if (centre.eq.'R') lattice = 3

          ! lattice is number of RL points per unit cell
          ! need to make this 3 if centre is 'R'
          do k = 1,lattice
            xff = x
            yff = y
            zff = z
            if ((k.eq.2).and.(centre.eq.'A')) then
              yff = y + 0.5d0
              zff = z + 0.5d0
            endif
            if ((k.eq.2).and.(centre.eq.'B')) then
              xff = x + 0.5d0
              zff = z + 0.5d0
            endif
            if ((k.eq.2).and.(centre.eq.'C')) then
              xff = x + 0.5d0
              yff = y + 0.5d0
            endif
            if ((k.eq.2).and.(centre.eq.'I')) then
              xff = x + 0.5d0
              yff = y + 0.5d0
              zff = z + 0.5d0
            endif
            if ((k.eq.2).and.(centre.eq.'F')) then
              xff = x + 0.5d0
              yff = y + 0.5d0
            endif
            if ((k.eq.3).and.(centre.eq.'F')) then
              xff = x + 0.5d0
              zff = z + 0.5d0
            endif
            if ((k.eq.4).and.(centre.eq.'F')) then
              yff = y + 0.5d0
              zff = z + 0.5d0
            endif

            if ((k.eq.2).and.(centre.eq.'R')) then
              xff = x + 2.0d0/3.0d0
              yff = y + 1.0d0/3.0d0
              zff = z + 1.0d0/3.0d0
            endif
            if ((k.eq.3).and.(centre.eq.'R')) then
              xff = x + 1.0d0/3.0d0
              yff = y + 2.0d0/3.0d0
              zff = z + 2.0d0/3.0d0
            endif

            if (ladjust_xyz_signs) then
              do while (xff>=1.0d0)
                xff = xff - 1.0d0
              end do
              do while (yff>=1.0d0)
                yff = yff - 1.0d0
              end do
              do while (zff>=1.0d0)
                zff = zff - 1.0d0
              end do
              do while (xff<0.0d0)
                xff = xff + 1.0d0
              end do
              do while (yff<0.0d0)
                yff = yff + 1.0d0
              end do
              do while (zff<0.0d0)
                zff = zff + 1.0d0
              end do
            end if
            okay = .true.
            if (natoms.gt.0) then
              do ij = 1,natoms
                dxf = xff - xf(ij) + 1.5d0
                dxf = dxf - aint(dxf) - 0.5d0
                dyf = yff - yf(ij) + 1.5d0
                dyf = dyf - aint(dyf) - 0.5d0
                dzf = zff - zf(ij) + 1.5d0
                dzf = dzf - aint(dzf) - 0.5d0
                dxo = cell(1,1)*dxf + cell(1,2)*dyf + cell(1,3)*dzf
                dyo = cell(2,1)*dxf + cell(2,2)*dyf + cell(2,3)*dzf
                dzo = cell(3,1)*dxf + cell(3,2)*dyf + cell(3,3)*dzf
                dr = dsqrt(dxo**2 + dyo**2 + dzo**2)
                if (dr<rsame) okay = .false.
              end do
            endif

            if (okay) then
              natoms = natoms + 1
              xf(natoms) = xff
              yf(natoms) = yff
              zf(natoms) = zff
              xo(natoms) = cell(1,1)*xff + cell(1,2)*yff + cell(1,3)*zff
              yo(natoms) = cell(2,1)*xff + cell(2,2)*yff + cell(2,3)*zff
              zo(natoms) = cell(3,1)*xff + cell(3,2)*yff + cell(3,3)*zff
              atom_name(natoms) = can(1:2)
!              if (luselabels) then
!                read(can(3:12),'(a)') catom_label(natoms)
!              else
!                read(can(3:12),*) atom_label(natoms)
!              end if
              ierror = 0
              read(can(3:8),*,iostat=ierror) atom_label(natoms)
              if (ierror/=0) atom_label(natoms) = 0
              reference_cell(natoms,:) = 0
              reference_number(natoms) = natoms
            endif
          end do  ! loop over k
        end do  ! loop over j
      end do ! loop over i
    end do ! loop over i2

    call assign_elements

! Count number of labels    
!   if (luselabels) then
! something here
!    end if

!    if (.not.luselabels) then
!      do i = 1,natoms
!        iel = atom_type(i)
!        n_labels(iel) = max(n_labels(iel),atom_label(i))
!      end do
!    end if

    do i = 1,natoms
      iel = atom_type(i)
      n_labels(iel) = max(n_labels(iel),atom_label(i))
    end do

    density = dble(natoms)/volume

    if (ldiag) then
      write(main_output,*) ' Number of sites: ',natoms
      write(main_output,*)
      write(main_output,'(a)') ' === Configuration in fractional coordinates ==='
      write(main_output,*) a,b,c
      write(main_output,*) alpha,beta,gamma
      write(main_output,'(a,i0)') 'Number of atoms = ',natoms
      do i = 1,natoms
!        if (luselabels) then
!          write(main_output,'(i0,2x,a2,2x,i0,2x,a10,2x,3f10.4)') i,atom_name(i),atom_type(i), &
!                catom_label(i),xf(i),yf(i),zf(i)
!        else
!          write(main_output,'(i0,2x,a2,2(2x,i0),2x,3f10.4)') i,atom_name(i),atom_type(i), &
!                atom_label(i),xf(i),yf(i),zf(i)
!        end if
        write(main_output,'(i0,2x,a2,2(2x,i0),2x,3f10.4)') i,atom_name(i),atom_type(i), &
                atom_label(i),xf(i),yf(i),zf(i)
      end do
      write(main_output,*)
      write(main_output,'(a)') ' === Configuration in orthogonal coordinates ==='
      write(main_output,*) a,b,c
      write(main_output,*) alpha,beta,gamma
      write(main_output,'(a,i0)') 'Number of atoms = ',natoms
      do i = 1,natoms
!        if (luselabels) then
!          write(main_output,'(i0,2x,a2,2x,i0,2x,a10,2x,3f10.4)') i,atom_name(i),atom_type(i), &
!               catom_label(i),xo(i),yo(i),zo(i)
!        else
!          write(main_output,'(i0,2x,a2,2(2x,i0),2x,3f10.4)') i,atom_name(i),atom_type(i), &
!                atom_label(i),xo(i),yo(i),zo(i)
!        end if
        write(main_output,'(i0,2x,a2,2(2x,i0),2x,3f10.4)') i,atom_name(i),atom_type(i), &
              atom_label(i),xo(i),yo(i),zo(i)
      end do
      write(main_output,*)
      do iel = -1,nelements
        if(n_elements(iel)>0) write(main_output,*) element(iel),n_elements(iel),n_labels(iel)
      end do
    end if

    return

  end subroutine generate_structure


!===============================================================================

  subroutine obtain_structure_cfg
! ===============================

!--------------------------------------------------------------------------
!
!     This subroutine reads in the configuration in a version 3 RMC
!     configuration (.cfg) file
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use annotations
    use channel_numbers

   implicit none

    integer :: i,j,k,n,iatom
    double precision :: xyz(3),pi,volume
    character(len=80) :: buffer

    pi = 2.0d0*dacos(0.0d0)

    read(ic,'(a)') buffer
    n = index(buffer,'(Version 3 format configuration file)')
    if (n==0) then
    write(6,'(a)') 'Configuration file does''t conform to version 3 format'
    write(6,'(a)') '... at least in the first line. Program ending.'
    stop
    end if

    read(ic,'(a)') buffer
    metadata_title = adjustl(trim(buffer))

    buffer=''
    do while (index(buffer,'moves generated, tried, accepted')==0)
    read(ic,'(a)') buffer
    end do
    read(buffer,*) moves_generated,moves_tried,moves_accepted

    buffer=''
    do while (index(buffer,'configurations saved')==0)
    read(ic,'(a)') buffer
    end do
    read(buffer,*) prior_saves

    buffer=''
    do while (index(buffer,'molecules (of all types)')==0)
    read(ic,'(a)') buffer
    end do
    read(buffer,*) natoms

    buffer=''
    do while (index(buffer,'types of molecule')==0)
    read(ic,'(a)') buffer
    end do
    read(buffer,*) ntypes

    allocate(numoftype(ntypes))
    allocate(aname(ntypes))

    buffer=''
    do while (index(buffer,'Defining vectors are:')==0)
    read(ic,'(a)') buffer
    end do
    do j = 1,3
    read(ic,*) cell(:,j)
    end do

    do i = 1,3
    do j = 1,3
      cell(i,j) = 2.0d0*cell(i,j)  ! Convert from cfg's half cell vectors
    end do
    end do

    a = dsqrt(cell(1,1)**2 + cell(2,1)**2 + cell(3,1)**2)
    b = dsqrt(cell(1,2)**2 + cell(2,2)**2 + cell(3,2)**2)
    c = dsqrt(cell(1,3)**2 + cell(2,3)**2 + cell(3,3)**2)
    alpha = (180.0d0/pi)*dacos((cell(1,2)*cell(1,3) + cell(2,2)*cell(2,3) + &
                cell(3,2)*cell(3,3))/(b*c))
    beta  = (180.0d0/pi)*dacos((cell(1,1)*cell(1,3) + cell(2,1)*cell(2,3) + &
                cell(3,1)*cell(3,3))/(a*c))
    gamma = (180.0d0/pi)*dacos((cell(1,1)*cell(1,2) + cell(2,1)*cell(2,2) + &
                cell(3,1)*cell(3,2))/(a*b))
    volume = cell(1,1)*cell(2,2)*cell(3,3) + cell(2,1)*cell(3,2)*cell(1,3) + &
           cell(3,1)*cell(1,2)*cell(2,3) - cell(1,1)*cell(3,2)*cell(2,3) - &
           cell(2,1)*cell(1,2)*cell(3,3) - cell(3,1)*cell(2,2)*cell(1,3)
!    write(6,*) volume
    natoms = 0
    do i = 1,ntypes
    buffer=''
    do while (index(buffer,'molecules of type')==0)
      read(ic,'(a)') buffer
    end do
    read(buffer,*) numoftype(i)
    natoms = natoms + numoftype(i)
    write(6,'(a,i0,a,i0,a)', advance="no") 'Please give name of atom type ',i,' (',numoftype(i),' atoms) : '
    read(5,*) aname(i)
    aname(i) = adjustl(aname(i))
    end do

    allocate(xf(natoms),yf(natoms),zf(natoms))
    allocate(xo(natoms),yo(natoms),zo(natoms))
    allocate(atom_name(natoms))
    allocate(atom_type(natoms))
    allocate(atom_label(natoms))
    allocate(reference_cell(natoms,3),reference_number(natoms))

    density = dble(natoms)/volume

    do i = 1,3
    read(ic,*)
    end do

    iatom = 0
    n_elements = 0
    do i = 1,ntypes
    do j = 1,numoftype(i)
      read(ic,*) xyz
      do k = 1,3
        xyz(k) = 0.5d0*(xyz(k)+1.0d0)
      end do
      iatom = iatom + 1
      atom_name(iatom) = aname(i)
      xf(iatom) = xyz(1) ; yf(iatom) = xyz(2) ; zf(iatom) = xyz(3)
      xo(iatom) = cell(1,1)*xyz(1) + cell(1,2)*xyz(2) + cell(1,3)*xyz(3)
      yo(iatom) = cell(2,1)*xyz(1) + cell(2,2)*xyz(2) + cell(2,3)*xyz(3)
      zo(iatom) = cell(3,1)*xyz(1) + cell(3,2)*xyz(2) + cell(3,3)*xyz(3)
      reference_cell(iatom,:) = 0
      reference_number(iatom) = iatom
    end do
    end do

    call assign_elements

    if (ldiag) then
      write(main_output,*)
      write(main_output,*) 'Output from obtain_struture_cfg'
      write(main_output,*) ' Number of sites: ',natoms
      write(main_output,*)
      write(main_output,*) a,b,c
      write(main_output,*) alpha,beta,gamma
      write(main_output,'(a)') ' === Original configuration in fractional coordinates ==='
      do i = 1,natoms
        write(main_output,'(i4,2x,a4,2x,i3,3f10.4)') i,atom_name(i),atom_type(i), &
                xf(i),yf(i),zf(i)
      end do
      write(main_output,*)
      write(main_output,*) a,b,c
      write(main_output,*) alpha,beta,gamma
      write(main_output,'(a)') ' === Original configuration in orthogonal coordinates ==='
      do i = 1,natoms
        write(main_output,'(i4,2x,a4,2x,i3,3f10.4)') i,atom_name(i),atom_type(i), &
               xo(i),yo(i),zo(i)
      end do
      write(main_output,*)
    end if

    return

  end subroutine obtain_structure_cfg

!===============================================================================

  subroutine get_cfg_spins
! ========================

!--------------------------------------------------------------------------
!
!     This subroutine reads in the magnetic spin configuration from a
!     version 3 format (.cfg) spin file
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use channel_numbers

   implicit none

    integer :: i

    open(imag,file=trim(filemag),status='old')
    read(imag,*)
    read(imag,*) nspins

    allocate(spin(nspins,3))

    do i = 1,nspins
    read(imag,*) spin(i,:)
    end do

    close(imag)

    return

  end subroutine get_cfg_spins

!===============================================================================

  subroutine obtain_structure_his
! ===============================

!--------------------------------------------------------------------------
!
!     This subroutine reads in the configuration in a version 3 RMC
!     histogram (.his) file
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use annotations
    use histogram_commons
    use channel_numbers

   implicit none

    integer :: i,j,k,n,jatoms,iatom,ios,ip
    double precision :: xyz(3),pi,volume
    character(len=80) :: buffer
    character(len=5) :: apairt,atemp1,atemp2

    pi = 2.0d0*dacos(0.0d0)

    read(ic,'(a)') buffer
    n = index(buffer,'RMCA (v3) intermediate (histogram) file')
    if (n==0) then
    write(6,'(a)') 'Histogram file does''t conform to version 3 format'
    write(6,'(a)') '... at least in the first line. Program ending.'
    stop
    end if

    read(ic,'(a)') buffer
    metadata_title = adjustl(trim(buffer)) ; config_title = metadata_title ; lm_title = .true.

    read(ic,*,IOSTAT=ios) moves_generated, moves_tried, moves_accepted, prior_saves
    if (ios>0) call read_error('moves_generated, moves_tried, moves_accepted, '// &
                              'prior_saves','his')
    if (ios<0) call end_of_file_error('moves_generated, moves_tried, moves_accepted, '// &
                                     'prior_saves','his')

    do j = 1,3
    read(ic,*,IOSTAT=ios) cell(:,j)
    if (ios>0) call read_error('cell','his')
    if (ios<0) call end_of_file_error('cell','his')
    end do

    do i = 1,3
    do j = 1,3
      cell(i,j) = 2.0d0*cell(i,j)  ! Convert from cfg's half cell vectors
    end do
    end do

    a = dsqrt(cell(1,1)**2 + cell(2,1)**2 + cell(3,1)**2)
    b = dsqrt(cell(1,2)**2 + cell(2,2)**2 + cell(3,2)**2)
    c = dsqrt(cell(1,3)**2 + cell(2,3)**2 + cell(3,3)**2)
    alpha = (180.0d0/pi)*dacos((cell(1,2)*cell(1,3) + cell(2,2)*cell(2,3) + &
                cell(3,2)*cell(3,3))/(b*c))
    beta  = (180.0d0/pi)*dacos((cell(1,1)*cell(1,3) + cell(2,1)*cell(2,3) + &
                cell(3,1)*cell(3,3))/(a*c))
    gamma = (180.0d0/pi)*dacos((cell(1,1)*cell(1,2) + cell(2,1)*cell(2,2) + &
                cell(3,1)*cell(3,2))/(a*b))
    volume = cell(1,1)*cell(2,2)*cell(3,3) + cell(2,1)*cell(3,2)*cell(1,3) + &
           cell(3,1)*cell(1,2)*cell(2,3) - cell(1,1)*cell(3,2)*cell(2,3) - &
           cell(2,1)*cell(1,2)*cell(3,3) - cell(3,1)*cell(2,2)*cell(1,3)

    read(ic,*) density,natoms,ntypes
    if (ios>0) call read_error('density,natoms,ntypes','his')
    if (ios<0) call end_of_file_error('density,natoms,ntypes','his')
    if (nint(density/volume)/=natoms) then
    density = dble(natoms)/volume
    write(buffer,*) density ; buffer = adjustl(buffer)
    write(6,'(a)') 'There is an inconsistency with density value in the .his file'
    write(6,'(a)') 'We are using a recomputed value of the density equal to '//trim(buffer)
    end if

    allocate(xf(natoms),yf(natoms),zf(natoms))
    allocate(xo(natoms),yo(natoms),zo(natoms))
    allocate(atom_name(natoms))
    allocate(atom_type(natoms))
    allocate(atom_label(natoms))

    allocate(numoftype(ntypes))
    allocate(aname(ntypes))

    read(ic,*,IOSTAT=ios) numoftype
    if (ios>0) call read_error('numoftype','his')
    if (ios<0) call end_of_file_error('numoftype','his')

    jatoms = 0
    do i = 1,ntypes
    jatoms = jatoms + numoftype(i)
    write(6,'(a,i0,a,i0,a)', advance="no") 'Please give name of atom type ',i,' (',numoftype(i),' atoms) : '
    read(5,*) aname(i)
    aname(i) = adjustl(aname(i))
    end do
    if (jatoms/=natoms) then
    write(6,'(a)') 'The configuration file doesn''t add up:'
    write(6,'(a)') 'Total number of atoms /= sum of numbers of each type'
    write(6,'(a)') '... so program has to end'
    stop
    end if

    iatom = 0
    n_elements = 0
    do i = 1,ntypes
    do j = 1,numoftype(i)
      read(ic,*,IOSTAT=ios) xyz
      if (ios>0) call read_error('xyz','his')
      if (ios<0) call end_of_file_error('xyz','his')
      do k = 1,3
        xyz(k) = 0.5d0*(xyz(k)+1.0d0)
      end do
      iatom = iatom + 1
      atom_name(iatom) = aname(i)
      xf(iatom) = xyz(1) ; yf(iatom) = xyz(2) ; zf(iatom) = xyz(3)
      xo(iatom) = cell(1,1)*xyz(1) + cell(1,2)*xyz(2) + cell(1,3)*xyz(3)
      yo(iatom) = cell(2,1)*xyz(1) + cell(2,2)*xyz(2) + cell(2,3)*xyz(3)
      zo(iatom) = cell(3,1)*xyz(1) + cell(3,2)*xyz(2) + cell(3,3)*xyz(3)
    end do
    end do

    call assign_elements

! Here we read the pair distribution functions
    read (ic,*,IOSTAT=ios) nr,dr
    if (ios>0) call read_error('nr,dr','his')
    if (ios<0) call end_of_file_error('nr,dr','his')
    npar = ntypes*(ntypes+1)/2
    allocate(histogram(nr,npar))
    do i = 1, npar
    read (ic,*,IOSTAT=ios) (histogram(j,i),j=1,nr)
    if (ios>0) call read_error('histogram','his')
    if (ios<0) call end_of_file_error('histogram','his')
    end do

! For the purpose of subsequent writing, construct the character array containing the names
! of the pairs of atoms in each PDF
    allocate(apairs(npar))
    ip = 0
    do i = 1,ntypes
    atemp1 = adjustl(aname(i))
    apairt(1:3) = atemp1(1:2)//'–'
    do j = i,ntypes
      atemp2 = adjustl(aname(j))
      apairt(4:5) = atemp2(1:2)
      if (apairt(2:2)==' ') apairt(2:) = apairt(3:)//' '
      ip = ip + 1
      apairs(ip) = apairt
    end do
    end do

! Finally we read
    read (ic,*,IOSTAT=ios) ncoord_0
    if (ncoord_0>0) then
    allocate(typec_0(ncoord_0))
    allocate(typen_0(ncoord_0))
    allocate(rcoord_0(2,ncoord_0))
    allocate(nneigh(natoms,ncoord_0))
    do i = 1, ncoord_0
      read(ic,*,IOSTAT=ios) typec_0(i),typen_0(i),rcoord_0(1,i), &
          rcoord_0(2,i),(nneigh(j,i),j=1,natoms)
      if (ios>0) call read_error('typec_0 etc','his')
      if (ios<0) call end_of_file_error('typec_0 etc','his')
    end do
    end if

  end subroutine obtain_structure_his

!===============================================================================

  subroutine obtain_structure_cssr
! ================================

!--------------------------------------------------------------------------
!
!     This subroutine reads in the configuration in a cssr file
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use channel_numbers

    implicit none

    integer :: i,n
    double precision :: xyz(3),volume,al,be,ga,s
    character(len=256) :: buffer

    read(ic,*) a,b,c
    read(ic,*) alpha,beta,gamma
    read(ic,*) natoms
!    read(ic,*)
!   Form the lattice vectors
    al = alpha/57.295779512d0
    be = beta/57.295779512d0
    ga = gamma/57.295779512d0
    s = 0.5*(al+be+ga)
!    volume = 2.0*a*b*c*sqrt(sin(s)*sin(s-al)*sin(s-be)*sin(s-ga))
    volume = a*b*c*sqrt(1.0d0-cos(2.0d0*al)-cos(2.0d0*be)-cos(2.0d0*ga)) + 2.0d0*(cos(al)*cos(be)*cos(ga))
    cell(1,1) = volume/(sin(al)*b*c)
    cell(2,1) = a*cos(ga)*sin(al)
    cell(3,1) = a*cos(be)
    cell(1,2) = 0.0d0
    cell(2,2) = b*sin(al)
    cell(3,2) = b*cos(al)
    cell(1,3) = 0.0d0
    cell(2,3) = 0.0d0
    cell(3,3) = c
    density = dble(natoms)/volume

    allocate(xf(natoms),yf(natoms),zf(natoms))
    allocate(xo(natoms),yo(natoms),zo(natoms))
    allocate(atom_name(natoms))
    allocate(atom_type(natoms))
    allocate(atom_label(natoms))
    allocate(reference_cell(natoms,3),reference_number(natoms))

    do i = 1,natoms
    read(ic,*) n,atom_name(i),xyz
    xf(i) = xyz(1) ; yf(i) = xyz(2) ; zf(i) = xyz(3)
    xo(i) = cell(1,1)*xyz(1) + cell(1,2)*xyz(2) + cell(1,3)*xyz(3)
    yo(i) = cell(2,1)*xyz(1) + cell(2,2)*xyz(2) + cell(2,3)*xyz(3)
    zo(i) = cell(3,1)*xyz(1) + cell(3,2)*xyz(2) + cell(3,3)*xyz(3)
    reference_cell(i,:) = 0
    reference_number(i) = i
    end do

    call assign_elements

    if (ldiag) then
    write(main_output,*)
    write(main_output,*) 'Output from obtain_struture_cssr'
    write(main_output,*) ' Number of sites: ',natoms
    write(main_output,*)
    write(main_output,'(a)') ' === Original configuration in fractional coordinates ==='
      write(main_output,*) a,b,c
      write(main_output,*) alpha,beta,gamma
    do i = 1,natoms
      write(main_output,'(i4,2x,a4,2x,i3,3f10.4)') i,atom_name(i),atom_type(i), &
              xf(i),yf(i),zf(i)
    end do
    write(main_output,*)
    write(main_output,'(a)') ' === Original configuration in orthogonal coordinates ==='
      write(main_output,*) a,b,c
      write(main_output,*) alpha,beta,gamma
    do i = 1,natoms
      write(main_output,'(i4,2x,a4,2x,i3,3f10.4)') i,atom_name(i),atom_type(i), &
             xo(i),yo(i),zo(i)
    end do
    write(main_output,*)
    end if

    return

  end subroutine obtain_structure_cssr

!===============================================================================

  subroutine obtain_structure_cell
! ================================

!--------------------------------------------------------------------------
!
!     This subroutine reads in the configuration in a Discus .cell file
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use channel_numbers

    implicit none

    integer :: i,n,j,ios,Reason, ncellx, ncelly, ncellz
    double precision :: xyz(3),volume,al,be,ga,s
    character(len=256) :: buffer, temp

    read(ic,'(a)') buffer
    read(ic,'(a)') buffer
    read(ic,'(a)') buffer
    j=index(buffer,',')
    do while (j>0)
        buffer = trim(buffer(1:j-1))//trim(buffer(j+1:))
        j=index(buffer,',')
    enddo
    
    read(buffer,*) temp, a,b,c,alpha,beta,gamma
    do while (trim(buffer).ne.'atoms')
        read(ic,'(a)') buffer
    enddo
!   Form the lattice vectors
    al = alpha/57.295779512d0
    be = beta/57.295779512d0
    ga = gamma/57.295779512d0
    s = 0.5*(al+be+ga)
!    volume = 2.0*a*b*c*sqrt(sin(s)*sin(s-al)*sin(s-be)*sin(s-ga))
    volume = a*b*c*sqrt(1.0d0-cos(2.0d0*al)-cos(2.0d0*be)-cos(2.0d0*ga)) + 2.0d0*(cos(al)*cos(be)*cos(ga))
    cell(1,1) = volume/(sin(al)*b*c)
    cell(2,1) = a*cos(ga)*sin(al)
    cell(3,1) = a*cos(be)
    cell(1,2) = 0.0d0
    cell(2,2) = b*sin(al)
    cell(3,2) = b*cos(al)
    cell(1,3) = 0.0d0
    cell(2,3) = 0.0d0
    cell(3,3) = c
    ! Count number of atoms in the file
    Reason = 0
    do while (Reason == 0)
        if (trim(buffer).ne.'') natoms = natoms + 1
       read(ic,*,IOSTAT=Reason) buffer
    enddo
    natoms = natoms - 1
!    close(ic)
    
    rewind(ic)
!    open(ic,file=trim(filename),form='formatted',status='old',iostat=ios)
    buffer = ''
    do while (trim(buffer).ne.'atoms')
        read(ic,'(a)') buffer
    enddo
    
    density = dble(natoms)/volume

    allocate(xf(natoms),yf(natoms),zf(natoms))
    allocate(xo(natoms),yo(natoms),zo(natoms))
    allocate(atom_name(natoms))
    allocate(atom_type(natoms))
    allocate(atom_label(natoms))
    allocate(reference_cell(natoms,3))
    allocate(reference_number(natoms))

    do i = 1,natoms
        read(ic,'(a)') buffer
        j=index(buffer,',')
        do while (j>0)
            buffer = trim(buffer(1:j-1))//trim(buffer(j+1:))
            j=index(buffer,',')
        enddo
        
    read(buffer,*) atom_name(i),xyz
    xf(i) = xyz(1) ; yf(i) = xyz(2) ; zf(i) = xyz(3)
    reference_cell(i,:) = 0
    reference_number(i) = i
    atom_label(i)       = 1
    end do
    
    ! Check supercell size
    ncellx = ceiling(maxval(xf))-floor(minval(xf))
    ncelly = ceiling(maxval(yf))-floor(minval(yf))
    ncellz = ceiling(maxval(zf))-floor(minval(zf))
    
    !   Form the lattice vectors
    a = a * ncellx
    b = b * ncelly
    c = c * ncellz
 !   volume = 2.0*a*b*c*sqrt(sin(s)*sin(s-al)*sin(s-be)*sin(s-ga))
    volume = a*b*c*sqrt(1.0d0-cos(2.0d0*al)-cos(2.0d0*be)-cos(2.0d0*ga)) + 2.0d0*(cos(al)*cos(be)*cos(ga))
    cell(1,1) = volume/(sin(al)*b*c)
    cell(2,1) = a*cos(ga)*sin(al)
    cell(3,1) = a*cos(be)
    cell(1,2) = 0.0d0
    cell(2,2) = b*sin(al)
    cell(3,2) = b*cos(al)
    cell(1,3) = 0.0d0
    cell(2,3) = 0.0d0
    cell(3,3) = c
    
    if (ncellx.ne.0.or.ncelly.ne.0.or.ncellz.ne.0) then
        write(6,'(a,i4,x,i4,x,i4)') 'Supercell size recognised from Discus is ', ncellx, ncelly, ncellz
        write(6,*) 'Check if this informaction is correct,'
        write(6,*) 'if not transfor the structure into a single unit cell'
    endif
    
    do i=1,natoms
        xf(i) = (xf(i) + floor(real(ncellx) / 2)) / real(ncellx)
        yf(i) = (yf(i) + floor(real(ncelly) / 2)) / real(ncelly)
        zf(i) = (zf(i) + floor(real(ncellz) / 2)) / real(ncellz)
        xo(i) = cell(1,1)*xf(i) + cell(1,2)*yf(i) + cell(1,3)*zf(i)
        yo(i) = cell(2,1)*xf(i) + cell(2,2)*yf(i) + cell(2,3)*zf(i)
        zo(i) = cell(3,1)*xf(i) + cell(3,2)*yf(i) + cell(3,3)*zf(i)
    enddo

    call assign_elements

    if (ldiag) then
    write(main_output,*)
    write(main_output,*) 'Output from obtain_struture_cell'
    write(main_output,*) ' Number of sites: ',natoms
    write(main_output,*)
    write(main_output,'(a)') ' === Original configuration in fractional coordinates ==='
    do i = 1,natoms
      write(main_output,'(i4,2x,a4,2x,i3,3f10.4)') i,atom_name(i),atom_type(i), &
              xf(i),yf(i),zf(i)
    end do
    write(main_output,*)
    write(main_output,'(a)') ' === Original configuration in orthogonal coordinates ==='
    do i = 1,natoms
      write(main_output,'(i4,2x,a4,2x,i3,3f10.4)') i,atom_name(i),atom_type(i), &
             xo(i),yo(i),zo(i)
    end do
    write(main_output,*)
    end if

    return

  end subroutine obtain_structure_cell

!===============================================================================  

  subroutine obtain_structure_rmc6f
! =================================

!--------------------------------------------------------------------------
!
!     This subroutine reads in the configuration in an rmc6f file
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use annotations
    use channel_numbers

    implicit none

    double precision :: s, volume, al,be,ga
    integer :: ios,i,n,itemp(4),ntemp
    logical :: lcell,lvectors
    character(len=132) :: buffer
    character(len=80) :: catom_types,cnumber_types

    ios = 0
    natoms = 0
    lcell = .false. ; lvectors = .false.
!    lsupercell = .false.
    metadata_title = ''
    metadata_owner = ''
    metadata_date = ''
    metadata_material = ''
    metadata_comment = ''
    metadata_source = ''
    read_loop: do while (ios==0)
      read(ic,'(a)',iostat=ios) buffer
      buffer = adjustl(buffer)
      n = index(buffer,':')
!!!      if (n==0) n = len(trim(buffer)) + 1
      if (n>0) then
        do i = 1,n-1                          ! this bit converts to uppercase
          if ((ichar(buffer(i:i))>=ichar('a')).and.(ichar(buffer(i:i))<=ichar('z'))) &
            buffer(i:i) = char(ichar(buffer(i:i)) - ichar('a') + ichar('A'))
        end do
      end if
      if (index(buffer,'NUMBER OF TYPES OF ATOMS')>0) then
        n = index(buffer,':') + 1
        read(buffer(n:),*) ntypes
        cycle read_loop
      end if
      if (index(buffer,'ATOM TYPES PRESENT')>0) then
        n = index(buffer,':') + 1
        read(buffer(n:),'(a)') catom_types
        cycle read_loop
      end if
      if (index(buffer,'NUMBER OF EACH ATOM TYPE')>0) then
        n = index(buffer,':') + 1
        read(buffer(n:),'(a)') cnumber_types
        cycle read_loop
      end if
      if (index(buffer,'METADATA TILE')>0) then
        n = index(buffer,':') + 1
        metadata_title = trim(adjustl(buffer(n:)))
        cycle read_loop
      end if
      if (index(buffer,'METADATA OWNER')>0) then
        n = index(buffer,':') + 1
        metadata_owner = trim(adjustl(buffer(n:)))
        cycle read_loop
      end if
      if (index(buffer,'METADATA DATE')>0) then
        n = index(buffer,':') + 1
        metadata_date = trim(adjustl(buffer(n:)))
        cycle read_loop
      end if
      if (index(buffer,'METADATA MATERIAL')>0) then
        n = index(buffer,':') + 1
        metadata_material = trim(adjustl(buffer(n:)))
        cycle read_loop
      end if
      if (index(buffer,'METADATA COMMENT')>0) then
        n = index(buffer,':') + 1
        metadata_comment = trim(adjustl(buffer(n:)))
        cycle read_loop
      end if
      if (index(buffer,'METADATA SOURCE')>0) then
        n = index(buffer,':') + 1
        metadata_source = trim(adjustl(buffer(n:)))
        cycle read_loop
      end if
      if (index(buffer,'NUMBER OF MOVES GENERATED')>0) then
        n = index(buffer,':') + 1
        read(buffer(n:),*) moves_generated
        cycle read_loop
      end if
      if (index(buffer,'NUMBER OF MOVES TRIED')>0) then
        n = index(buffer,':') + 1
        read(buffer(n:),*) moves_tried
        cycle read_loop
      end if
      if (index(buffer,'NUMBER OF MOVES ACCEPTED')>0) then
        n = index(buffer,':') + 1
        read(buffer(n:),*) moves_accepted
        cycle read_loop
      end if
      if (index(buffer,'NUMBER OF PRIOR CONFIGURATION SAVES')>0) then
        n = index(buffer,':') + 1
        read(buffer(n:),*) prior_saves
        cycle read_loop
      end if
      if (index(buffer,'NUMBER OF ATOMS')>0) then
        n = index(buffer,':') + 1
        read(buffer(n:),*) natoms
        cycle read_loop
      end if
!      if (index(buffer,'NUMBER OF SPINS')>0) then
!        n = index(buffer,':') + 1
!        read(buffer(n:),*) bragg_nspins
!        cycle read_loop
!      end if
      if (index(buffer,'SUPERCELL DIMENSIONS')>0) then
        n = index(buffer,':') + 1
        ! >>>>>>>>>>>> Yuanpeng -> when the input file 
        ! is '.rmc6f' file, the supercell dimension 
        ! will be read into an independent array. Since 
        ! otherwise the 'ncell' read in from here will 
        ! overwrite the one read in from the arguments.
        read(buffer(n:),*) ncell_rmc6f
        ! <<<<<<<<<<<< Yuanpeng finishes here <<<<<<<<<
        cycle read_loop
      end if
      if (index(buffer,'CELL (ANG/DEG)')>0) then
        n = index(buffer,':') + 1
        read(buffer(n:),*) a,b,c,alpha,beta,gamma
        lcell = .true.
        cycle read_loop
      end if
      if (index(buffer,'LATTICE VECTORS (ANG)')>0) then
        read(ic,*,iostat=ios) cell
        if (ios==0) lvectors = .true.
        cycle read_loop
      end if
      if_atoms: if (index(buffer,'ATOMS')>0) then
! Note that we are reading lines that look like
!     1   Na  [1]  0.000000    0.000000    0.000000     1   0   0   0
! The first integer is optional, so we need to test for its existence.
! The integer in square brackets is also optional, and we test for the bracket.
! Moreover, the last four integers are also optional, and again we need to check they exist
! There are two blocks of text, one for the data file defining a value for natoms and one
! where we need to count the number of atoms. We use the technique of looking for leading white
! space and then removing it using the ADJUSTL function.
        if_natoms: if (natoms>0) then
          allocate(xf(natoms),yf(natoms),zf(natoms))
          allocate(xo(natoms),yo(natoms),zo(natoms))
          allocate(atom_name(natoms))
          allocate(atom_type(natoms))
          allocate(atom_label(natoms))
!!!          allocate(atom_number(natoms))
!!!          allocate(reference_cell(natoms,3),reference_number(natoms))
          atom_label = 0   ! assign for the case where no label is given
          atom_loop_1: do i = 1,natoms
            read(ic,'(a)',iostat=ios) buffer
            if (ios/=0) then
              write(6,'(a)') 'Exit> End of rmc6f file encountered during atom read'
              stop
            end if
            buffer = adjustl(buffer)                  ! first check for atom number
            if ( (ichar(buffer(1:1))>=ichar('0')).and.(ichar(buffer(1:1))<=ichar('9')) ) then
!!!              read(buffer,*) atom_number(i)           ! read atom number
              n = index(buffer,' ')                   ! move past the atom number
              buffer = adjustl(buffer(n:))
!!!            else
!!!              atom_number(i) = i                      ! otherwise assign atom number
            end if
            atom_name(i) = buffer(1:2)           ! read element label
            buffer = buffer(3:)                       ! move past the element symbol
!***
!***  Define element_symbol(:), atom_number(:), atom_label(:)
!***
            buffer = adjustl(buffer)
            if (buffer(1:1)=='[') then                ! this bit to read atom label
              n = index(buffer,']')
              read(buffer(2:n-1),*) atom_label(i)
              buffer = buffer(n+1:)
            end if
            read(buffer,*) xf(i),yf(i),zf(i)          ! read x,y,z
            buffer = adjustl(buffer)
            n = index(buffer,' ') ; buffer = adjustl(buffer(n:))   ! move to the point on
            n = index(buffer,' ') ; buffer = adjustl(buffer(n:))   ! the line that contains
            n = index(buffer,' ') ; buffer = adjustl(buffer(n:))   ! integers
            if (len(trim(buffer))==0) then   ! case where there are no integers
              ntemp = 0
              itemp = 0
            else
              read(buffer,*) itemp(1)
              n = index(buffer,' ') ; buffer = adjustl(buffer(n:))
              if (len(trim(buffer))==0) then ! case where there is 1 integer
                ntemp = 1
              else
                read(buffer,*) itemp(2)
                n = index(buffer,' ') ; buffer = adjustl(buffer(n:))
                read(buffer,*) itemp(3)
                n = index(buffer,' ') ; buffer = adjustl(buffer(n:))
                if (len(trim(buffer))==0) then   ! case where there are 3 integers
                  ntemp = 3
                else                             ! case where there are 4 integers
                  read(buffer,*) itemp(4)
                  ntemp = 4
                end if
              end if
            end if
            if (i==1) then
              ! >>>>>>>>>>>> Yuanpeng -> debugging >>>>>>>>>>>>>>>>>>>>>>>>>>>
              ! The condition 'ntemp==0' was added to deal with the >>>>>>>>>>
              ! situation where no reference cell information is given >>>>>>>
              ! in the '.rmc6f' input file. Previously when 'ntemp==0', >>>>>>
              ! 'reference_number' and 'reference_cell' was not allocated >>>>
              ! but later on the 'expand_structure' subroutine will try to >>>
              ! get access to the 'reference_number' and 'reference_cell' >>>>
              ! arrays which then causes the segmentation fault. >>>>>>>>>>>>>
              if ((ntemp==1).or.(ntemp==4).or.(ntemp==0)) allocate (reference_number(natoms))
              if ((ntemp==4).or.(ntemp==0)) allocate (reference_cell(natoms,3))
              ! <<<<<<<<<<<<< Yuanpeng's editing finishes here <<<<<<<<<<<<<<<
            end if
            if ((ntemp==1).or.(ntemp==4)) reference_number(i) = itemp(1)
            ! >>>>>>>>>>>>>>>> Yuanpeng -> debugging >>>>>>>>>>>>>>>>>>>>>>>>>
            ! The following statement for assigning values for >>>>>>>>>>>>>>>
            ! 'reference_number' and 'reference_cell' arrays was copied >>>>>>
            ! from the subroutine like 'obtain_structure_cell'. >>>>>>>>>>>>>>
            if (ntemp==0) then
              reference_cell(i,:) = 0
              reference_number(i) = i
            end if
            ! <<<<<<<<<<<<<<<< Yuanpeng's editing finishes here <<<<<<<<<<<<<<
            if (ntemp==3) reference_cell(i,:) = itemp(1:3)
            if (ntemp==4) reference_cell(i,:) = itemp(2:4)
          end do atom_loop_1
        else
          allocate(xf(1),yf(1),zf(1))
          allocate(atom_name(1))
!!!          allocate(atom_number(1))
          allocate(atom_type(1))
          natoms = 0
          atom_loop_2: do while (ios==0)
            read(ic,'(a)',iostat=ios) buffer
            if (ios/=0) exit atom_loop_2
            natoms = natoms + 1 ; i = natoms
            if (i>1) call extend_arrays(i)
            buffer = adjustl(buffer)                  ! first check for atom number
            if ((ichar(buffer(1:1))>=ichar('0')).and.(ichar(buffer(1:1))<=ichar('9'))) then
!!!              read(buffer,*) atom_number(i)
              n = index(buffer,' ')                   ! move past the atom number
              buffer = adjustl(buffer(n:))
!!!            else
!!!              atom_number(i) = i                      ! otherwise assign atom number
            end if
            atom_name(i) = buffer(1:2)
            buffer = buffer(3:)                       ! move past the element label
            buffer = adjustl(buffer)
            if (buffer(1:1)=='[') then                ! this bit to read atom label
              n = index(buffer,']')
!!              read(buffer(1:n-1),*) atom_type(i)
              buffer = buffer(n+1:)
            end if
            read(buffer,*) xf(i),yf(i),zf(i)
            buffer = adjustl(buffer)
            n = index(buffer,' ') ; buffer = adjustl(buffer(n:))   ! move to the point on
            n = index(buffer,' ') ; buffer = adjustl(buffer(n:))   ! the line that contains
            n = index(buffer,' ') ; buffer = adjustl(buffer(n:))   ! integers
            if (len(trim(buffer))==0) then   ! case where there are no integers
              ntemp = 0
              itemp = 0
            else
              read(buffer,*) itemp(1)
              n = index(buffer,' ') ; buffer = adjustl(buffer(n:))
              if (len(trim(buffer))==0) then ! case where there is 1 integer
                ntemp = 1
              else
                read(buffer,*) itemp(2)
                n = index(buffer,' ') ; buffer = adjustl(buffer(n:))
                read(buffer,*) itemp(3)
                n = index(buffer,' ') ; buffer = adjustl(buffer(n:))
                if (len(trim(buffer))==0) then   ! case where there are 3 integers
                  ntemp = 3
                else                             ! case where there are 4 integers
                  read(buffer,*) itemp(4)
                  ntemp = 4
                end if
              end if
            end if
            if ((ntemp==1).or.(ntemp==4)) reference_number(i) = itemp(1)
            if (ntemp==3) reference_cell(i,:) = itemp(1:3)
            if (ntemp==4) reference_cell(i,:) = itemp(2:4)
          end do atom_loop_2
        end if if_natoms
      end if if_atoms
     end do read_loop

    call assign_elements

! check that we were given information about the lattice
   if ((.not.lcell).and.(.not.lvectors)) then
     write(6,'(a)') 'info> no lattice information given in rmc6f file'
     stop
   end if
! if lattice vectors are not provided, obtain from the lattice parameters
   if (.not.lvectors) then
     al = alpha*dasin(1.0d0)/90.0d0
     be = beta*dasin(1.0d0)/90.0d0
     ga = gamma*dasin(1.0d0)/90.0d0
     s = 0.5d0*(alpha+beta+gamma)
!     volume = 2.0d0*a*b*c*dsqrt(dsin(s)*dsin(s-alpha)*dsin(s-beta)*dsin(s-gamma))
     volume = a*b*c*sqrt(1.0d0-cos(2.0d0*al)-cos(2.0d0*be)-cos(2.0d0*ga)) + 2.0d0*(cos(al)*cos(be)*cos(ga))
     cell(1,1) = volume/dsin(alpha)/b/c
     cell(2,1) = a*dcos(gamma)*dsin(alpha)
     cell(3,1) = a*dcos(beta)
     cell(1,2) = 0.0d0
     cell(2,2) = b*dsin(alpha)
     cell(3,2) = b*dcos(alpha)
     cell(1,3) = 0.0d0
     cell(2,3) = 0.0d0
     cell(3,3) = c
!     alpha = alpha*90.0d0/dasin(1.0d0)
!     beta = beta*90.0d0/dasin(1.0d0)
!     gamma = gamma*90.0d0/dasin(1.0d0)
   else
     volume = cell(1,1)*cell(2,2)*cell(3,3) + cell(1,2)*cell(2,3)*cell(3,1) + &
              cell(1,3)*cell(2,1)*cell(3,2) - cell(1,1)*cell(2,3)*cell(3,2) - &
              cell(1,2)*cell(2,1)*cell(3,3) - cell(1,3)*cell(2,2)*cell(3,1)
   end if

   do i = 1,natoms
      xo(i) = cell(1,1)*xf(i) + cell(1,2)*yf(i) + cell(1,3)*zf(i)
      yo(i) = cell(2,1)*xf(i) + cell(2,2)*yf(i) + cell(2,3)*zf(i)
      zo(i) = cell(3,1)*xf(i) + cell(3,2)*yf(i) + cell(3,3)*zf(i)
   end do

   density = dble(natoms)/volume

   close(ic)

    if (ldiag) then
      write(main_output,*)
      write(main_output,*) 'Output from obtain_struture_rmc6f'
      write(main_output,*) ' Number of sites: ',natoms
      write(main_output,*)
      write(main_output,'(a)') ' === Original configuration in fractional coordinates ==='
      do i = 1,natoms
        write(main_output,'(i4,2x,a4,2x,i3,3f10.4)') i,atom_name(i),atom_type(i), &
                xf(i),yf(i),zf(i)
      end do
      write(main_output,*)
      write(main_output,'(a)') ' === Original configuration in orthogonal coordinates ==='
      do i = 1,natoms
        write(main_output,'(i4,2x,a4,2x,i3,3f10.4)') i,atom_name(i),atom_type(i), &
               xo(i),yo(i),zo(i)
      end do
      write(main_output,*)
    end if


    return

  end subroutine obtain_structure_rmc6f

!===============================================================================

  subroutine obtain_structure_dlpoly
! ==================================

!--------------------------------------------------------------------------
!
!     This subroutine reads in the configuration from a DL_POLY CONFIG of
!     REVCON file
!
!--------------------------------------------------------------------------

   use structure_data
   use arguments
   use channel_numbers

   implicit none

   integer :: i,ilines,ireadflag
   double precision :: volume,x,y,z,rvolume,rcell(3,3)
   double precision :: calpha, cbeta, cgamma
   character(len=80) :: cbuffer,title

   read(ic,'(a)') title
   read(ic,*) ilines    ! Note, only read the first integer from this line
   do i=1,3 ; read(ic,*) cell(i,:) ; end do
   volume = cell(1,1)*cell(2,2)*cell(3,3) + cell(1,2)*cell(2,3)*cell(3,1) + &
            cell(1,3)*cell(2,1)*cell(3,2) - cell(1,1)*cell(2,3)*cell(3,2) - &
            cell(1,3)*cell(2,2)*cell(3,1) - cell(1,2)*cell(2,1)*cell(3,3)
   a = sqrt(cell(1,1)**2 + cell(1,2)**2 + cell(1,3)**2)
   b = sqrt(cell(2,1)**2 + cell(2,2)**2 + cell(2,3)**2)
   c = sqrt(cell(3,1)**2 + cell(3,2)**2 + cell(3,3)**2)
   calpha = (cell(2,1)*cell(3,1) + cell(2,2)*cell(3,2) + cell(2,3)*cell(3,3))/(b*c)
   cbeta  = (cell(1,1)*cell(3,1) + cell(1,2)*cell(3,2) + cell(1,3)*cell(3,3))/(a*c)
   cgamma = (cell(1,1)*cell(2,1) + cell(1,2)*cell(2,2) + cell(1,3)*cell(2,3))/(a*b)
   alpha = acos(calpha)*90.0d0/asin(1.0d0)
   beta  = acos(cbeta)*90.0d0/asin(1.0d0)
   gamma = acos(cgamma)*90.0d0/asin(1.0d0)

   volume = cell(1,1)*cell(2,2)*cell(3,3) + cell(1,2)*cell(2,3)*cell(3,1) + &
            cell(1,3)*cell(2,1)*cell(3,2) - cell(1,1)*cell(2,3)*cell(3,2) - &
            cell(1,3)*cell(2,2)*cell(3,1) - cell(1,2)*cell(2,1)*cell(3,3)

   rcell(1,1) = cell(2,2)*cell(3,3) - cell(2,3)*cell(3,2) ! correct
   rcell(1,2) = cell(2,3)*cell(3,1) - cell(2,1)*cell(3,3) ! correct
   rcell(1,3) = cell(2,1)*cell(3,2) - cell(2,2)*cell(3,1) ! correct
   rcell(2,1) = cell(3,2)*cell(1,3) - cell(3,3)*cell(1,2) ! correct
   rcell(2,2) = cell(3,3)*cell(1,1) - cell(3,1)*cell(1,3) ! correct
   rcell(2,3) = cell(3,1)*cell(1,2) - cell(3,2)*cell(1,1) ! correct
   rcell(3,1) = cell(1,2)*cell(2,3) - cell(1,3)*cell(2,2) ! correct
   rcell(3,2) = cell(1,3)*cell(2,1) - cell(1,1)*cell(2,3) ! correct
   rcell(3,3) = cell(1,1)*cell(2,2) - cell(1,2)*cell(2,1) ! correct

   rcell = rcell/volume

   rvolume = rcell(1,1)*rcell(2,2)*rcell(3,3) + rcell(1,2)*rcell(2,3)*rcell(3,1) + &
             rcell(1,3)*rcell(2,1)*rcell(3,2) - rcell(1,1)*rcell(2,3)*rcell(3,2) - &
             rcell(1,3)*rcell(2,2)*rcell(3,1) - rcell(1,2)*rcell(2,1)*rcell(3,3)

   natoms = 0
   ireadflag = 0
   readloop: do while (ireadflag==0)
     read(ic,'(a)',iostat=ireadflag) cbuffer
     if (ireadflag/=0) exit readloop
     natoms = natoms + 1
     read(ic,*) x,y,z
     if (ilines>=1) read(ic,*)
     if (ilines==2) read(ic,*)
   end do readloop
   rewind(ic)

   allocate(xf(natoms),yf(natoms),zf(natoms))
   allocate(xo(natoms),yo(natoms),zo(natoms))
   allocate(atom_name(natoms))
   allocate(atom_type(natoms))
!   if (luselabels) then
!     allocate(catom_label(natoms))
!   else
!     allocate(atom_label(natoms))
!   end if
   allocate(atom_label(natoms))
   allocate(reference_cell(natoms,3),reference_number(natoms))

   density = dble(natoms)/volume

   do i = 1,5 ; read(ic,*) ; end do
   do i = 1,natoms
     read(ic,'(a)') cbuffer
     cbuffer = adjustl(cbuffer)
     atom_name(i) = cbuffer(1:2)
     read(ic,*) x,y,z
     xo(i) = x ; yo(i) = y ; zo(i) = z
     xf(i) = rcell(1,1)*x + rcell(1,2)*y + rcell(1,3)*z
     yf(i) = rcell(2,1)*x + rcell(2,2)*y + rcell(2,3)*z
     zf(i) = rcell(3,1)*x + rcell(3,2)*y + rcell(3,3)*z
     reference_cell(i,:) = 0
     reference_number(i) = i
     if (ilines>=1) read(ic,*)
     if (ilines==2) read(ic,*)
   end do

   call assign_elements

   if (ldiag) then
   write(main_output,*)
   write(main_output,*) 'Output from obtain_struture_dlpoly'
   write(main_output,*) ' Number of atoms: ',natoms
   write(main_output,*)
   write(main_output,'(a)') ' === Original configuration in fractional coordinates ==='
   do i = 1,natoms
     write(main_output,'(i4,2x,a4,2x,i3,3f10.4)') i,atom_name(i),atom_type(i),xf(i),yf(i),zf(i)
   end do
   write(main_output,*)
   write(main_output,'(a)') ' === Original configuration in orthogonal coordinates ==='
   do i = 1,natoms
     write(main_output,'(i4,2x,a4,2x,i3,3f10.4)') i,atom_name(i),atom_type(i),xo(i),yo(i),zo(i)
   end do
   write(main_output,*)
   end if

   return

  end subroutine obtain_structure_dlpoly


!===============================================================================


  subroutine obtain_structure_atomeye
! ===================================

!--------------------------------------------------------------------------
!
!     This subroutine reads in the configuration from an atomeye file
!
!--------------------------------------------------------------------------

   use structure_data
   use arguments
   use channel_numbers

   implicit none

   integer :: i,ilines,ireadflag,n,nentry
   double precision :: volume,xyz(3),rvolume,rcell(3,3)
   double precision :: calpha, cbeta, cgamma, scale_factor
   character(len=80) :: cbuffer,title

   do
     read(ic,'(a)') cbuffer
     if (len_trim(cbuffer)==0) cycle
     cbuffer = adjustl(cbuffer)
     n = index(cbuffer,'=')
     select case(trim(cbuffer(1:n-1)))
       case('Number of particles')
         read(cbuffer(n+1:),*) natoms
       case('A')
         read(cbuffer(n+1:),*) scale_factor
       case('H0(1,1)')
         read(cbuffer(n+1:),*) cell(1,1)         
       case('H0(2,1)')
         read(cbuffer(n+1:),*) cell(2,1)         
       case('H0(3,1)')
         read(cbuffer(n+1:),*) cell(3,1)         
       case('H0(1,2)')
         read(cbuffer(n+1:),*) cell(1,2)         
       case('H0(2,2)')
         read(cbuffer(n+1:),*) cell(2,2)         
       case('H0(3,2)')
         read(cbuffer(n+1:),*) cell(3,2)         
       case('H0(1,3)')
         read(cbuffer(n+1:),*) cell(1,3)         
       case('H0(2,3)')
         read(cbuffer(n+1:),*) cell(2,3)         
       case('H0(3,3)')
         read(cbuffer(n+1:),*) cell(3,3)     
       case('.NO_VELOCITY.')
         continue
       case('entry_count')
         read(cbuffer(n+1:),*) nentry
         exit
     end select
   end do  
   if (nentry>3) then 
     do i = 4,nentry
       read(ic,*)
     end do
   end if
   
   cell = cell*scale_factor

   volume = cell(1,1)*cell(2,2)*cell(3,3) + cell(1,2)*cell(2,3)*cell(3,1) + &
            cell(1,3)*cell(2,1)*cell(3,2) - cell(1,1)*cell(2,3)*cell(3,2) - &
            cell(1,3)*cell(2,2)*cell(3,1) - cell(1,2)*cell(2,1)*cell(3,3)
   a = sqrt(cell(1,1)**2 + cell(1,2)**2 + cell(1,3)**2)
   b = sqrt(cell(2,1)**2 + cell(2,2)**2 + cell(2,3)**2)
   c = sqrt(cell(3,1)**2 + cell(3,2)**2 + cell(3,3)**2)
   calpha = (cell(2,1)*cell(3,1) + cell(2,2)*cell(3,2) + cell(2,3)*cell(3,3))/(b*c)
   cbeta  = (cell(1,1)*cell(3,1) + cell(1,2)*cell(3,2) + cell(1,3)*cell(3,3))/(a*c)
   cgamma = (cell(1,1)*cell(2,1) + cell(1,2)*cell(2,2) + cell(1,3)*cell(2,3))/(a*b)
   alpha = acos(calpha)*90.0d0/asin(1.0d0)
   beta  = acos(cbeta)*90.0d0/asin(1.0d0)
   gamma = acos(cgamma)*90.0d0/asin(1.0d0)

   rcell(1,1) = cell(2,2)*cell(3,3) - cell(2,3)*cell(3,2) ! correct
   rcell(1,2) = cell(2,3)*cell(3,1) - cell(2,1)*cell(3,3) ! correct
   rcell(1,3) = cell(2,1)*cell(3,2) - cell(2,2)*cell(3,1) ! correct
   rcell(2,1) = cell(3,2)*cell(1,3) - cell(3,3)*cell(1,2) ! correct
   rcell(2,2) = cell(3,3)*cell(1,1) - cell(3,1)*cell(1,3) ! correct
   rcell(2,3) = cell(3,1)*cell(1,2) - cell(3,2)*cell(1,1) ! correct
   rcell(3,1) = cell(1,2)*cell(2,3) - cell(1,3)*cell(2,2) ! correct
   rcell(3,2) = cell(1,3)*cell(2,1) - cell(1,1)*cell(2,3) ! correct
   rcell(3,3) = cell(1,1)*cell(2,2) - cell(1,2)*cell(2,1) ! correct

   rcell = rcell/volume

   rvolume = rcell(1,1)*rcell(2,2)*rcell(3,3) + rcell(1,2)*rcell(2,3)*rcell(3,1) + &
             rcell(1,3)*rcell(2,1)*rcell(3,2) - rcell(1,1)*rcell(2,3)*rcell(3,2) - &
             rcell(1,3)*rcell(2,2)*rcell(3,1) - rcell(1,2)*rcell(2,1)*rcell(3,3)

   allocate(xf(natoms),yf(natoms),zf(natoms))
   allocate(xo(natoms),yo(natoms),zo(natoms))
   allocate(atom_name(natoms))
   allocate(atom_type(natoms))
!   if (luselabels) then
!     allocate(catom_label(natoms))
!   else
!     allocate(atom_label(natoms))
!   end if
   allocate(atom_label(natoms))
   allocate(reference_cell(natoms,3),reference_number(natoms))

   density = dble(natoms)/volume
   
   do i = 1,natoms
     read(ic,*) 
     read(ic,'(a)') cbuffer
     cbuffer = adjustl(cbuffer)
     atom_name(i) = cbuffer(1:2)
     read(ic,*) xyz
    xf(i) = xyz(1) ; yf(i) = xyz(2) ; zf(i) = xyz(3)
    xo(i) = cell(1,1)*xyz(1) + cell(1,2)*xyz(2) + cell(1,3)*xyz(3)
    yo(i) = cell(2,1)*xyz(1) + cell(2,2)*xyz(2) + cell(2,3)*xyz(3)
    zo(i) = cell(3,1)*xyz(1) + cell(3,2)*xyz(2) + cell(3,3)*xyz(3)
!      xo(i) = x ; yo(i) = y ; zo(i) = z
!      xf(i) = rcell(1,1)*x + rcell(1,2)*y + rcell(1,3)*z
!      yf(i) = rcell(2,1)*x + rcell(2,2)*y + rcell(2,3)*z
!      zf(i) = rcell(3,1)*x + rcell(3,2)*y + rcell(3,3)*z
     reference_cell(i,:) = 0
     reference_number(i) = i
   end do

   call assign_elements

!    natoms = 0
!    ireadflag = 0
!    readloop: do while (ireadflag==0)
!      read(ic,'(a)',iostat=ireadflag) cbuffer
!      if (ireadflag/=0) exit readloop
!      natoms = natoms + 1
!      read(ic,*) x,y,z
!      if (ilines>=1) read(ic,*)
!      if (ilines==2) read(ic,*)
!    end do readloop
!    rewind(ic)

   if (ldiag) then
   write(main_output,*)
   write(main_output,*) 'Output from obtain_struture_atomeye'
   write(main_output,*) ' Number of atoms: ',natoms
   write(main_output,*)
   write(main_output,'(a)') ' === Original configuration in fractional coordinates ==='
   do i = 1,natoms
     write(main_output,'(i4,2x,a4,2x,i3,3f10.4)') i,atom_name(i),atom_type(i),xf(i),yf(i),zf(i)
   end do
   write(main_output,*)
   write(main_output,'(a)') ' === Original configuration in orthogonal coordinates ==='
   do i = 1,natoms
     write(main_output,'(i4,2x,a4,2x,i3,3f10.4)') i,atom_name(i),atom_type(i),xo(i),yo(i),zo(i)
   end do
   write(main_output,*)
   end if

   return

  end subroutine obtain_structure_atomeye


!===============================================================================


  subroutine obtain_structure_www
! ===============================

!--------------------------------------------------------------------------
!
!     Parses the structure information from a www file
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use channel_numbers

    implicit none

    integer :: iatom,ierror,n,i
    double precision :: box,xyz(3)
    character(len=80) :: buffer
    character(len=8) :: catom
    
    nsym = 1

! Count the number of atoms, checking that the file looks ok
    read(ic,'(a)') buffer
    n = index(buffer,'=') 
    if (n==0) then
      write(6,'(a)') 'The first line of the www file seems wrong'
      stop
    end if
    read(ic,*)
    natoms = 0  ;  ierror = 0
    checkfile: do
      read(ic,'(a)',iostat=ierror) buffer
      if (len_trim(buffer)==0) exit checkfile
      if (ierror/=0) exit checkfile
      read(buffer,*,iostat=ierror) iatom
      if (ierror/=0) then
        write(6,'(a)') 'Something is wrong with the format of the www file'
        stop
      end if
      if (iatom==natoms) then
        natoms = natoms + 1
        cycle checkfile
      end if
      if (iatom==0.and.natoms>0) then
        if (iatom>2*natoms) then
          write(6,'(a)') 'Something is wrong with the number of atoms in the www file'
          stop
        end if
      end if
    end do checkfile
    
    allocate(xf(natoms),yf(natoms),zf(natoms))
    allocate(xo(natoms),yo(natoms),zo(natoms))
    allocate(atom_name(natoms))
    allocate(atom_type(natoms))
    allocate(atom_label(natoms))
    allocate(www_neighbours(natoms))

    rewind(ic)
    read(ic,'(a)') buffer
    n = index(buffer,'=') 
    read(buffer(n+1:),*) box
    a = box*2.0d0
    b = a  ;  c = a
    alpha = 90.0d0  ; beta = 90.0d0  ;  gamma = 90.0d0
    cell = 0.0d0
    cell(1,1) = a  ;  cell(2,2) = a  ;  cell(3,3) = a
    
    read(ic,*)

    do i = 1,natoms
      read(ic,*) iatom,xyz
      xo(i) = xyz(1) + box
      yo(i) = xyz(2) + box
      zo(i) = xyz(3) + box
      xf(i) = xo(i)/a
      yf(i) = yo(i)/a
      zf(i) = zo(i)/a
    end do 
    
    atom_name = 'Si'

    call assign_elements

    do i = 1,natoms
      read(ic,'(a)') www_neighbours(i)
    end do 

    if (ldiag) then
      write(main_output,*)
      write(main_output,*) 'Output from obtain_struture_www'
      write(main_output,*) ' Number of atoms: ',natoms
      write(main_output,*)
      write(main_output,'(a)') ' === Original configuration in fractional coordinates ==='
      do i = 1,natoms
        write(main_output,'(i0,2x,a4,2x,i3,3f10.4)') i,atom_name(i),atom_type(i),xf(i),yf(i),zf(i)
      end do
      write(main_output,*)
      write(main_output,'(a)') ' === Original configuration in orthogonal coordinates ==='
      do i = 1,natoms
        write(main_output,'(i0,2x,a4,2x,i3,3f10.4)') i,atom_name(i),atom_type(i),xo(i),yo(i),zo(i)
      end do
      write(main_output,*)
    end if

    return  

  end subroutine obtain_structure_www


!===============================================================================


  subroutine extend_arrays(n)
! ===========================

!..Extends dimensions of atom arrays during reading from configurations. This accounts for
!  the case where the input configuration file doesn't contain the number of atoms and thus
!  the number has to be determined through a line-by-line count during the reading stage.

  use structure_data
  use arguments

  implicit none

  integer :: n,m,shape1(1),shape2(2),iatoms,iatom_type, &
             ireference_cell,ireference_number,iatom_name
  double precision, allocatable :: tatoms(:,:)
  integer, allocatable :: treference_cell(:,:),treference_number(:),tatom_type(:), &
                          tatom_label(:)
  character(len=10), allocatable :: ctatom_label(:)
  character(len=4), allocatable :: tatom_name(:)

  m = n-1

 !         allocate(xf(natoms),yf(natoms),zf(natoms))
 !         allocate(xo(natoms),yo(natoms),zo(natoms))
 !         allocate(atom_name(natoms))
 !         allocate(atom_type(natoms))


! check that the bookwork is in order
  shape1 = shape(xf)
  iatoms = shape1(1)
  shape1 = shape(atom_name)
  iatom_name = shape1(1)
  shape1 = shape(atom_type)
  iatom_type = shape1(1)
  if (iatoms/=m) then
    write(6,'(a)') 'Exit> extend_arrays has found problem with book-keeping of array atoms'
    stop
  end if
  if (iatom_name/=m) then
    write(6,'(a)') 'Exit> extend_arrays has found problem with book-keeping of array atom_label'
    stop
  end if
  allocate(tatoms(3,m)) ; tatoms(1,:) = xf ; tatoms(2,:) = yf ; tatoms(3,:) = zf
  allocate(tatom_name(m)) ; tatom_name = atom_name
  allocate(tatom_type(m)) ; tatom_type = atom_type
!  if (luselabels) then
!    allocate(ctatom_label(m)) ; ctatom_label = catom_label
!  else
!    allocate(tatom_label(m)) ; tatom_label = atom_label
!  end if
  allocate(tatom_label(m)) ; tatom_label = atom_label
  if (allocated(xf)) then
    deallocate(xf) ; allocate(xf(n)) ; xf(1:m) = tatoms(1,:)
  end if
  if (allocated(yf)) then
    deallocate(yf) ; allocate(yf(n)) ; yf(1:m) = tatoms(2,:)
  end if
  if (allocated(zf)) then
    deallocate(zf) ; allocate(zf(n)) ; zf(1:m) = tatoms(3,:)
  end if
  if (allocated(atom_name)) then
    deallocate(atom_name) ; allocate(atom_name(n)) ; atom_name(1:m) = tatom_name
  end if
  if (allocated(atom_type)) then
    deallocate(atom_type) ; allocate(atom_type(n)) ; atom_type(1:m) = tatom_type
  end if
  if (allocated(atom_label)) then
    deallocate(atom_label) ; allocate(atom_label(n)) ; atom_label(1:m) = tatom_label
  end if
  if (allocated(catom_label)) then
    deallocate(catom_label) ; allocate(catom_label(n)) ; catom_label(1:m) = ctatom_label
  end if
  if (allocated(tatoms)) deallocate(tatoms)
  if (allocated(tatom_name)) deallocate(tatom_name)
  if (allocated(tatom_type)) deallocate(tatom_type)
  if (allocated(tatom_label)) deallocate(tatom_label)
  if (allocated(ctatom_label)) deallocate(ctatom_label)
  if (allocated(reference_cell)) then
    shape2 = shape(reference_cell)
    ireference_cell = 1 + shape2(1)*shape2(2)/3
    if (ireference_cell/=m) then
      write(6,'(a)') &
      'Exit> extend_arrays has found problem with book-keeping of array reference_cell'
      stop
    end if
    allocate(treference_cell(3,m))
    if (allocated(reference_cell)) deallocate(reference_cell)
    allocate(reference_cell(3,n))
    reference_cell(:,1:m) = treference_cell
    if (allocated(treference_cell)) deallocate(treference_cell)
  end if
  if (allocated(reference_number)) then
    shape1 = shape(reference_number)
    ireference_number = shape1(1)
    if (ireference_number/=m) then
      write(6,'(a)') &
      'Exit> extend_arrays has found problem with book-keeping of array reference_number'
      stop
    end if
    allocate(treference_number(m))
    if (allocated(reference_number)) deallocate(reference_number)
    allocate(reference_number(n))
    reference_number(1:m) = treference_number
    if (allocated(treference_number)) deallocate(treference_number)
  end if

  return

  end subroutine extend_arrays

!===============================================================================

  subroutine obtain_structure_crystal
! ===================================

!--------------------------------------------------------------------------
!
!     Parses the structure information from a crystal format file
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use channel_numbers

    implicit none

    integer :: i,j,n,ios,jatoms,iatom
    double precision :: volume,pi
    double precision, allocatable :: xft(:),yft(:),zft(:)

    ios = 0

! Read multipliers
    read(ic,*) ncell

! Read lattice vectors
    do j = 1,3
    read(ic,*,IOSTAT=ios) cell(:,j)
    if (ios>0) call read_error('cell','his')
    if (ios<0) call end_of_file_error('cell','his')
    end do

    pi = 2.0d0*dacos(0.0d0)
    a = dsqrt(cell(1,1)**2 + cell(2,1)**2 + cell(3,1)**2)
    b = dsqrt(cell(1,2)**2 + cell(2,2)**2 + cell(3,2)**2)
    c = dsqrt(cell(1,3)**2 + cell(2,3)**2 + cell(3,3)**2)
    alpha = (180.0d0/pi)*dacos((cell(1,2)*cell(1,3) + cell(2,2)*cell(2,3) + &
                cell(3,2)*cell(3,3))/(b*c))
    beta  = (180.0d0/pi)*dacos((cell(1,1)*cell(1,3) + cell(2,1)*cell(2,3) + &
                cell(3,1)*cell(3,3))/(a*c))
    gamma = (180.0d0/pi)*dacos((cell(1,1)*cell(1,2) + cell(2,1)*cell(2,2) + &
                cell(3,1)*cell(3,2))/(a*b))
    volume = cell(1,1)*cell(2,2)*cell(3,3) + cell(2,1)*cell(3,2)*cell(1,3) + &
           cell(3,1)*cell(1,2)*cell(2,3) - cell(1,1)*cell(3,2)*cell(2,3) - &
           cell(2,1)*cell(1,2)*cell(3,3) - cell(3,1)*cell(2,2)*cell(1,3)

    read(ic,*) ntypes
    allocate(numoftype(ntypes))
    allocate(aname(ntypes))

    natoms = 0
    do i = 1,ntypes
    read(ic,*) numoftype(i)
    n = numoftype(i)
    if (natoms>0) then
      allocate(xft(natoms),yft(natoms),zft(natoms))
      xft = xf; yft = yf; zft = zf
      if (allocated(xf)) deallocate(xf,yf,zf)
    end if
    allocate(xf(n+natoms),yf(n+natoms),zf(n+natoms))
    if (natoms>0) then
      xf(1:natoms) = xft ; yf(1:natoms) = yft ; zf(1:natoms) = zft
    end if
    do j = 1,n
      read(ic,*) xf(natoms+j),yf(natoms+j),zf(natoms+j)
    end do
    if (allocated(xft)) deallocate(xft,yft,zft)
    natoms = natoms + n
    end do

    allocate(xo(natoms),yo(natoms),zo(natoms))
    do i = 1,natoms
    xo(i) = cell(1,1)*xf(i) + cell(1,2)*yf(i) + cell(1,3)*zf(i)
    yo(i) = cell(2,1)*xf(i) + cell(2,2)*yf(i) + cell(2,3)*zf(i)
    zo(i) = cell(3,1)*xf(i) + cell(3,2)*yf(i) + cell(3,3)*zf(i)
    end do

    allocate(atom_name(natoms))
    allocate(atom_type(natoms))
!    if (luselabels) then
!      allocate(atom_label(natoms))
!    else
!      allocate(catom_label(natoms))
!    end if
    allocate(atom_label(natoms))
    allocate(reference_cell(natoms,3),reference_number(natoms))

    jatoms = 0
    do i = 1,ntypes
    jatoms = jatoms + numoftype(i)
    write(6,'(a,i0,a,i0,a)', advance="no") 'Please give name of atom type ',i,' (',numoftype(i),' atoms) : '
    read(5,*) aname(i)
    aname(i) = adjustl(aname(i))
    end do
    if (jatoms/=natoms) then
    write(6,'(a)') 'The configuration file doesn''t add up:'
    write(6,'(a)') 'Total number of atoms /= sum of numbers of each type'
    write(6,'(a)') '... so program has to end'
    stop
    end if

    iatom = 0
    do i = 1,ntypes
    do j = 1,numoftype(i)
      iatom = iatom + 1
      atom_name(iatom) = aname(i)
      reference_cell(iatom,:) = 0
      reference_number(iatom) = iatom
    end do
    end do

    density = dble(natoms)/volume

    call assign_elements

    if (ldiag) then
    write(main_output,*)
    write(main_output,*) 'Output from obtain_struture_crystal'
    write(main_output,*) ' Number of sites: ',natoms
    write(main_output,*)
    write(main_output,'(a)') ' === Original configuration in fractional coordinates ==='
    do i = 1,natoms
      write(main_output,'(i4,2x,a4,2x,i3,3f10.4)') i,atom_name(i),atom_type(i), &
              xf(i),yf(i),zf(i)
    end do
    write(main_output,*)
    write(main_output,'(a)') ' === Original configuration in orthogonal coordinates ==='
    do i = 1,natoms
      write(main_output,'(i4,2x,a4,2x,i3,3f10.4)') i,atom_name(i),atom_type(i), &
             xo(i),yo(i),zo(i)
    end do
    write(main_output,*)
    end if

    return

   end subroutine obtain_structure_crystal

!===============================================================================

  subroutine sort_structure
! =========================

!--------------------------------------------------------------------------
!
!     This subroutine sorts the list of atoms by element type.
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use channel_numbers

    implicit none

    integer :: i,j,iatom,iel,ilabel
    integer, allocatable :: attype(:),atlabel(:),rnumt(:),rcellt(:,:)
    double precision, allocatable :: xft(:),yft(:),zft(:),xot(:),yot(:),zot(:)
    character (len=4), allocatable :: atname(:)
    character (len=10), allocatable :: catlabel(:)

    allocate(attype(natoms),atlabel(natoms),catlabel(natoms),rnumt(natoms),rcellt(natoms,3))
    allocate(xft(natoms),yft(natoms),zft(natoms))
    allocate(xot(natoms),yot(natoms),zot(natoms))
    allocate(atname(natoms))

    xft = xf(1:natoms) ; yft = yf(1:natoms) ; zft = zf(1:natoms)
    xot = xo(1:natoms) ; yot = yo(1:natoms) ; zot = zo(1:natoms)
    attype = atom_type(1:natoms) ; atname = atom_name(1:natoms)
!    if (luselabels) then
!      catlabel = catom_label(1:natoms)
!    else
!      atlabel = atom_label(1:natoms)
!    end if
    atlabel = atom_label(1:natoms)
    rnumt = reference_number(1:natoms)
    rcellt = reference_cell(1:natoms,:)
    iatom = 0
    if (ldiag) write(main_output,'(/a)') 'Output from sort_structure'
    if (ldiag) write(main_output,'(a/)') '=========================='
    if (ldiag) write(main_output,'(a,i0)') 'Number of atoms at start = ',natoms
    do iel = -1,nelements
      if ((iel==0).and.ldelete) cycle
      if (n_elements(iel)==0) cycle
      if (n_labels(iel)==0) then
        do i = 1,natoms
          if (attype(i)==iel) then
            iatom = iatom + 1
            xf(iatom) = xft(i) ; yf(iatom) = yft(i) ; zf(iatom) = zft(i)
            xo(iatom) = xot(i) ; yo(iatom) = yot(i) ; zo(iatom) = zot(i)
            atom_name(iatom) = atname(i) ; atom_type(iatom) = attype(i)
!            if (luselabels) then
!              catom_label(iatom) = catlabel(i)
!            else
!              atom_label(iatom) = atlabel(i)
!            end if
            atom_label(iatom) = atlabel(i)
            reference_number(iatom) = rnumt(i)
            reference_cell(iatom,:) = rcellt(i,:)
          end if
        end do
      else
        do ilabel = 1,n_labels(iel)
          do i = 1,natoms
            if ((attype(i)==iel).and.(atlabel(i)==ilabel)) then
              iatom = iatom + 1
              xf(iatom) = xft(i) ; yf(iatom) = yft(i) ; zf(iatom) = zft(i)
              xo(iatom) = xot(i) ; yo(iatom) = yot(i) ; zo(iatom) = zot(i)
              atom_name(iatom) = atname(i) ; atom_type(iatom) = attype(i)
!              if (luselabels) then
!                catom_label(iatom) = catlabel(i)
!              else
!                atom_label(iatom) = atlabel(i)
!              end if
              atom_label(iatom) = atlabel(i)
              reference_number(iatom) = rnumt(i)
              reference_cell(iatom,:) = rcellt(i,:)
            end if
          end do
        end do      
      end if
      if (ldiag) write(main_output,*) 'iel,iatom',iel,iatom
    end do

    if (ldiag) then
      write(main_output,*)
      write(main_output,*) 'Old, new atom counts: ',natoms,iatom
      write(main_output,*)
    end if

    natoms = iatom
    if (ldelete) n_elements(0) = 0

    if (ldiag) then
      write(main_output,'(a)') ' === Sorted configuration in fractional coordinates ==='
      do i = 1,natoms
!        if (luselabels) then
!          write(main_output,'(i4,2x,a4,2x,i0,2x,a10,3f10.4,4i4)') i,atom_name(i),atom_type(i), &
!                catom_label(i),xf(i),yf(i),zf(i),reference_number(i),reference_cell(i,:)
!        else
!          write(main_output,'(i4,2x,a4,2x,i0,2x,i0,3f10.4,4i4)') i,atom_name(i),atom_type(i), &
!                atom_label(i),xf(i),yf(i),zf(i),reference_number(i),reference_cell(i,:)
!        end if
        write(main_output,'(i4,2x,a4,2x,i0,2x,i0,3f10.4,4i4)') i,atom_name(i),atom_type(i), &
              atom_label(i),xf(i),yf(i),zf(i),reference_number(i),reference_cell(i,:)
      end do
      write(main_output,*)
      write(main_output,'(a)') ' === Sorted configuration in orthogonal coordinates ==='
      do i = 1,natoms
        write(main_output,'(i4,2x,a4,2x,i3,3f10.4,4i4)') i,atom_name(i),atom_type(i), &
               xo(i),yo(i),zo(i),reference_number(i),reference_cell(i,:)
      end do
      write(main_output,*)
      do i = -1,nelements
        if (n_elements(i)==0) cycle
        write(main_output,'(a,i0,a,i0)') 'The configuration contains ',n_elements(i), &
                                ' atoms of element type ',i
      end do
      write(main_output,*)
    end if

    if (allocated(attype)) deallocate(attype)
    if (allocated(atlabel)) deallocate(atlabel)
    if (allocated(rnumt)) deallocate(rnumt)
    if (allocated(rcellt)) deallocate(rcellt)
    if (allocated(xft)) deallocate(xft,yft,zft)
    if (allocated(xot)) deallocate(xot,yot,zot)
    if (allocated(atname)) deallocate(atname)

    return

  end subroutine sort_structure

!===============================================================================

  subroutine order_structure
! ==========================

!--------------------------------------------------------------------------
!
!     This subroutine sorts the list of atoms by user selected order.
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use channel_numbers

    implicit none

    logical :: lnotrecoginsed
    integer :: i,iatom,iel,iell,itype,n,nchecked,ilabel
    integer, allocatable :: attype(:),atlabel(:),rnumt(:),rcellt(:,:),nlabel(:)
    double precision, allocatable :: xft(:),yft(:),zft(:),xot(:),yot(:),zot(:)
    character (len=10), allocatable :: catlabel(:)
    character (len=4), allocatable :: atname(:)
    character (len=2) :: catoml
    character (len=20) :: catomU
    logical :: lassigned
    logical, allocatable :: lused(:)

    allocate(attype(natoms),atlabel(natoms),catlabel(natoms),rnumt(natoms),rcellt(natoms,3))
    allocate(xft(natoms),yft(natoms),zft(natoms))
    allocate(xot(natoms),yot(natoms),zot(natoms))
    allocate(atname(natoms))
    allocate(lused(-1:nelements))

    lused = .false.
    lnotrecoginsed = .false.

    xft = xf(1:natoms) ; yft = yf(1:natoms) ; zft = zf(1:natoms)
    xot = xo(1:natoms) ; yot = yo(1:natoms) ; zot = zo(1:natoms)
    attype = atom_type(1:natoms) ; atname = atom_name(1:natoms)
!    if (luselabels) then
!      catlabel = catom_label(1:natoms)
!    else
!      atlabel = atom_label(1:natoms)
!    end if
    atlabel = atom_label(1:natoms)
    rnumt = reference_number(1:natoms)
    rcellt = reference_cell(1:natoms,:)
    iatom = 0
    iell = 0

! . Counts the number of element types being used, if this information was not previously
! . available
    if (ntypes==0) then
      do i = -1,nelements
        if (n_elements(i)>0) ntypes = ntypes + 1
      end do
    end if
    allocate(norder(ntypes)) ; norder = 0

! . Determines which atoms are present in the configuration if this information was not
! . previously available
    if (.not.allocated(numoftype)) then
      allocate(numoftype(ntypes))
      numoftype = 0
      i = 0
      do iel = -1,nelements
        if (n_elements(iel)>0) then
          i = i + 1
          numoftype(i) = n_elements(iel)
        end if
      end do
    end if

! . There are two options here. First is to receive the atom list via the -order [..]
! . flag. The second is for when the [...] isn't given, in which case the order is
! . requested at the terminal.

    if(lorder_list) then

! ... Here we extract the atoms to be swapped from the input list, and if any
! ... don't work, we simply extract the remainder of the list from the atoms not
! ... requested for ordering.
      itype = 0
      nchecked = 0
      corder_list = adjustl(corder_list)
      call remove_spaces(trim(corder_list))
      
      provided_atom_list_loop: do while (len_trim(corder_list)>0)
        read(corder_list(1:2),'(a)') catoml
        if (nchecked>ntypes) exit provided_atom_list_loop
        n = index(corder_list,' ')
        corder_list = adjustl(corder_list(n:))
! This bit checks that the atom names are real
        lassigned = .false.

        atom_check_loop: do while ((.not.lassigned).and.(.not.catoml==""))
            
! ..... Loop for elements with two letters in name
          sub_loopa: do iel = -1,nelements
            if ((catoml==element(iel)).or.(catoml==elementuc(iel)).or.(catoml==elementlc(iel))) then
              if ((n_elements(iel)>0).and.(.not.(lused(iel)))) then
                lassigned = .true.
                lused(iel) = .true.
                itype = itype + 1
                if (itype>ntypes) exit provided_atom_list_loop
                norder(itype) = iel
                nchecked = nchecked + 1
                exit atom_check_loop
              end if
            end if
          end do sub_loopa
          
          nchecked = nchecked + 1
          
          if (.not.lassigned) then
              lnotrecoginsed = .true.
              cycle provided_atom_list_loop ! cycle if atom is not specified    
          end if
       
        end do atom_check_loop
   
      end do provided_atom_list_loop

!     Exiting when one of the elements was not recoginsed or not present in the structure      
      if ((lnotrecoginsed).or.(nchecked/=ntypes)) call order_error

      sub_loopb: do iel = 1,nelements
        if ((n_elements(iel)>0).and.(.not.(lused(iel)))) then
          itype = itype + 1
          norder(itype) = iel
        end if
      end do sub_loopb

! Needs cheching this part of changes - Jun-23-2017
! >>>>>>>>>>>>>>>> Yuanpeng -> The following lines added by WS >>>>>>>>>>>>
! >>>>>>>>>>>>>>>> was commented out by Yuanpeng. The reason is >>>>>>>>>>>
! >>>>>>>>>>>>>>>> that the following codes will arrange atoms in >>>>>>>>>
! >>>>>>>>>>>>>>>> a way that cannot be dealt with properly by GASP! >>>>>>
! >>>>>>>>>>>>>>>> Currently not sure whether removing these lines >>>>>>>>
! >>>>>>>>>>>>>>>> will break any functionalities. But at least GASP >>>>>>
! >>>>>>>>>>>>>>>> is now working! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! WS_add here is the main loop for atom order
! I will add another level to the loop through the atlabel[]
! . Determines number of differetn atom sites for each element
!    allocate(nlabel(ntypes))
!    nlabel = 0
!    do itype = 1,ntypes
!        iel = norder(itype)
!        do i = 1,natoms
!          if (attype(i)==iel) then
!              nlabel(itype) = max(nlabel(itype),atlabel(i))        
!          endif
!        enddo
!    enddo
! <<<<<<<<<<<<<<<<<< Yuanpeng -> Editing finished here <<<<<<<<<<<<<<<<<<<

      numoftype = 0
      iatom = 0
      do itype = 1,ntypes
! >>>>>>>>>>>> Yuanpeng -> The following line commented out >>>>>>>>>>>>>>
! >>>>>>>>>>>> to account for the above mentioned bug with GASP >>>>>>>>>>
!       do ilabel = 1,nlabel(itype)
        iel = norder(itype)
        do i = 1,natoms
! >>>>>>>>>>>> Yuanpeng -> The following line commented out >>>>>>>>>>>>>>
! >>>>>>>>>>>> to account for the above mentioned bug with GASP >>>>>>>>>>
!          if (attype(i)==iel.and.atlabel(i)==ilabel) then
	   if (attype(i)==iel) then
            iatom = iatom + 1
            xf(iatom) = xft(i) ; yf(iatom) = yft(i) ; zf(iatom) = zft(i)
            xo(iatom) = xot(i) ; yo(iatom) = yot(i) ; zo(iatom) = zot(i)
            atom_name(iatom) = atname(i) ; atom_type(iatom) = attype(i)
!            if (luselabels) then
!              catom_label(iatom) = catlabel(i)
!            else
!              atom_label(iatom) = atlabel(i)
!            end if
            atom_label(iatom) = atlabel(i)
            reference_number(iatom) = rnumt(i)
            reference_cell(iatom,:) = rcellt(i,:)
            numoftype(itype) = numoftype(itype) + 1
          end if
        end do
! >>>>>>>>>>>> Yuanpeng -> The following line commented out >>>>>>>>>>>>>>
! >>>>>>>>>>>> to account for the above mentioned bug with GASP >>>>>>>>>>
!       end do
      end do
! <<<<<<<<<<<< Yuanpeng -> Editing finished here <<<<<<<<<<<<<<<<<<<<<<<<<

    else

! ... The following loop asks for the names of atoms.
! ... It does some error checking in order:
! ... a) Checks that the atom label is real, and asks for the information again if it decides
!       that the atom label isn't real and isn't in the configuration.
! ... b) Checks that the atom label hasn't been used more than once.
      itype = 0
      write(6,'(a)') 'Please give names of atoms in prefered order for the configuration file...'
      atom_list_loop: do while (itype<ntypes)
        itype = itype + 1
        write(6,'(a,i0,a)', advance="no") 'Name of atom for order position ',itype,': '
        read(5,'(a)') catoml
        catoml = adjustl(trim(catoml))

! This bit checks that the atom names are real
        lassigned = .false.
        element_loop: do while (.not.lassigned)
! ..... Loop for elements with two letters in name
          sub_loop1: do iel = -1,nelements
            if ((catoml==element(iel)).or.(catoml==elementuc(iel)).or.(catoml==elementlc(iel))) then
              if (n_elements(iel)>0) then
                lassigned = .true.
                iell = iel
                exit element_loop
              end if
            end if
          end do sub_loop1
! ....... Ask if element can't be assigned
          write(6,'(3a)') 'Element type for atom label ',trim(catoml),' not recognised in configuration...'
          write(6,'(a)', advance="no") 'Please give element symbol again: '
          read(5,'(a)') catoml
          catoml = adjustl(trim(catoml))
        end do element_loop

! ..... This bit checks that the atom hasn't been used, then if life is okay it puts the atoms
! ..... in the desired order
        if (lused(iell)) then
          write(6,'(a)') 'Detected duplicate atom, need to try again'
          cycle atom_list_loop
        else
          lused(iell) = .true.
          do i = 1,natoms
            if (attype(i)==iell) then
              iatom = iatom + 1
              xf(iatom) = xft(i) ; yf(iatom) = yft(i) ; zf(iatom) = zft(i)
              xo(iatom) = xot(i) ; yo(iatom) = yot(i) ; zo(iatom) = zot(i)
              atom_name(iatom) = atname(i) ; atom_type(iatom) = attype(i)
!              if (luselabels) then
!                catom_label(iatom) = catlabel(i)
!              else
!                atom_label(iatom) = atlabel(i)
!              end if
              atom_label(iatom) = atlabel(i)
              reference_number(iatom) = rnumt(i)
              reference_cell(iatom,:) = rcellt(i,:)
            end if
          end do
        end if
      end do atom_list_loop

    end if

    if (ldiag) then
    write(main_output,'(/a)') 'Output from order_structure'
    write(main_output,'(a/)') '==========================='
    write(main_output,*)
    write(main_output,*) 'Old, new atom counts: ',natoms,iatom
    write(main_output,*)
    write(main_output,'(a)') ' === Sorted configuration in fractional coordinates ==='
    do i = 1,natoms
      write(main_output,'(i4,2x,a4,2x,i3,3f10.4,4i4)') i,atom_name(i),atom_type(i), &
              xf(i),yf(i),zf(i),reference_number(i),reference_cell(i,:)
    end do
    write(main_output,*)
    write(main_output,'(a)') ' === Sorted configuration in orthogonal coordinates ==='
    do i = 1,natoms
      write(main_output,'(i4,2x,a4,2x,i3,3f10.4,4i4)') i,atom_name(i),atom_type(i), &
             xo(i),yo(i),zo(i),reference_number(i),reference_cell(i,:)
    end do
    write(main_output,*)
    end if

    if (allocated(attype)) deallocate(attype)
    if (allocated(atlabel)) deallocate(atlabel)
    if (allocated(rnumt)) deallocate(rnumt)
    if (allocated(rcellt)) deallocate(rcellt)
    if (allocated(xft)) deallocate(xft,yft,zft)
    if (allocated(xot)) deallocate(xot,yot,zot)
    if (allocated(atname)) deallocate(atname)
    if (allocated(nlabel)) deallocate(nlabel)

    return

  end subroutine order_structure


!===============================================================================


    subroutine shift_coordinates
!   ============================

    use structure_data

    implicit none

    xf = 2.0d0 + xf
    xf = xf - int(xf)
    yf = 2.0d0 + yf
    yf = yf - int(yf)
    zf = 2.0d0 + zf
    zf = zf - int(zf)
    
    xo = cell(1,1)*xf + cell(1,2)*yf + cell(1,3)*zf
    yo = cell(2,1)*xf + cell(2,2)*yf + cell(2,3)*zf
    zo = cell(3,1)*xf + cell(3,2)*yf + cell(3,3)*zf

   end subroutine shift_coordinates

    subroutine Uiso_structure
!   ============================

    use structure_data
    use random
    use utilities
    use arguments, only : Uiso_factor

    implicit none
    
    integer :: i
    double precision :: xUiso, yUiso, zUiso  ! x,y and z atom displacement in Cartezian coordinates
    double precision :: inverse_cell(3,3)
    
    inverse_cell = matrix33_inverse(cell)
	
    do i = 1, natoms

            xUiso = ZBQLNOR(0.0d0, Uiso(atom_type(i)) * Uiso_factor)
            yUiso = ZBQLNOR(0.0d0, Uiso(atom_type(i)) * Uiso_factor)
            zUiso = ZBQLNOR(0.0d0, Uiso(atom_type(i)) * Uiso_factor)
            xf(i) = xf(i) + xUiso * inverse_cell(1,1) + yUiso * inverse_cell(1,2) + zUiso * inverse_cell(1,3)
            yf(i) = yf(i) + xUiso * inverse_cell(2,1) + yUiso * inverse_cell(2,2) + zUiso * inverse_cell(2,3)
            zf(i) = zf(i) + xUiso * inverse_cell(3,1) + yUiso * inverse_cell(3,2) + zUiso * inverse_cell(3,3)           
    enddo
    

   end subroutine Uiso_structure


!===============================================================================


    subroutine expand_structure
!   ===========================

!--------------------------------------------------------------------------
!
!     This subroutine expands the list of atoms by adding extra lattice
!     vectors to the basis set.
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use channel_numbers

    implicit none

    integer :: i,j,k,ii,jj,kk,n,ia,newatoms,ncellt(3)
    integer, allocatable :: atom_typet(:),atom_labelt(:),rnumt(:),rcellt(:,:)
    double precision,allocatable :: xft(:),yft(:),zft(:),xot(:),yot(:),zot(:)
    character(len=10), allocatable :: catom_labelt(:)
    character(len=4), allocatable :: atom_namet(:)

    allocate(xft(natoms),yft(natoms),zft(natoms))
    allocate(xot(natoms),yot(natoms),zot(natoms))
    allocate(atom_namet(natoms))
    allocate(atom_typet(natoms))
!    if (luselabels) then
!      allocate(catom_labelt(natoms))
!    else
!      allocate(atom_labelt(natoms))
!    end if
    allocate(atom_labelt(natoms))
    allocate(rnumt(natoms),rcellt(natoms,3))
    xft = xf(1:natoms)
    yft = yf(1:natoms)
    zft = zf(1:natoms)
    xot = xo(1:natoms)
    yot = yo(1:natoms)
    zot = zo(1:natoms)
    atom_namet = atom_name(1:natoms)
    atom_typet = atom_type(1:natoms)
!    if (luselabels) then
!      catom_labelt = catom_label(1:natoms)
!    else
!      atom_labelt = atom_label(1:natoms)
!    end if
    atom_labelt = atom_label(1:natoms)
    rnumt = reference_number(1:natoms)
    rcellt = reference_cell(1:natoms,:)
    if (allocated(xf)) deallocate(xf,yf,zf)
    if (allocated(xo)) deallocate(xo,yo,zo)
    if (allocated(atom_name)) deallocate(atom_name)
    if (allocated(atom_type)) deallocate(atom_type)
    if (allocated(atom_label)) deallocate(atom_label)
    if (allocated(catom_label)) deallocate(catom_label)
    if (allocated(reference_number)) deallocate(reference_number)
    if (allocated(reference_cell)) deallocate(reference_cell)
    newatoms = natoms*ncell(1)*ncell(2)*ncell(3)
    do i = -1,nelements
      n_elements(i) = n_elements(i)*ncell(1)*ncell(2)*ncell(3)
    end do
    allocate(xf(newatoms),yf(newatoms),zf(newatoms))
    allocate(xo(newatoms),yo(newatoms),zo(newatoms))
    allocate(atom_name(newatoms))
    allocate(atom_type(newatoms))
!    if (luselabels) then
!      allocate(catom_label(newatoms))
!    else
!      allocate(atom_label(newatoms))
!    end if
    allocate(atom_label(newatoms))
    allocate(reference_number(newatoms))
    allocate(reference_cell(newatoms,3))
    xf(1:natoms) = xft
    yf(1:natoms) = yft
    zf(1:natoms) = zft
    xo(1:natoms) = xot
    yo(1:natoms) = yot
    zo(1:natoms) = zot
    atom_name(1:natoms) = atom_namet
    atom_type(1:natoms) = atom_typet
!    if (luselabels) then
!      catom_label(1:natoms) = catom_labelt
!    else
!      atom_label(1:natoms) = atom_labelt
!    end if
    atom_label(1:natoms) = atom_labelt
    reference_number(1:natoms) = rnumt
    reference_cell(1:natoms,:) = rcellt
    if (allocated(xft)) deallocate(xft,yft,zft)
    if (allocated(xot)) deallocate(xot,yot,zot)
    if (allocated(atom_namet)) deallocate(atom_namet)
    if (allocated(atom_typet)) deallocate(atom_typet)
    if (allocated(atom_labelt)) deallocate(atom_labelt)
    if (allocated(catom_labelt)) deallocate(catom_labelt)
    if (allocated(rnumt)) deallocate(rnumt)
    if (allocated(rcellt)) deallocate(rcellt)
    
    if (lreduce) then
      ncellt = ncell
      ncell = 1
    end if

    n = natoms
    do i = 1,ncell(1)
    ii = i-1
      do j = 1,ncell(2)
      jj = j-1
        loopk: do k = 1,ncell(3)
          if ((i.eq.1).and.(j.eq.1).and.(k.eq.1)) cycle loopk
          kk = k-1
          do ia = 1,natoms
            n = n + 1
            xf(n) = xf(ia) + dble(ii)
            yf(n) = yf(ia) + dble(jj)
            zf(n) = zf(ia) + dble(kk)
            xo(n) = xo(ia) + dble(ii)*cell(1,1) + dble(jj)*cell(1,2) + dble(kk)*cell(1,3)
            yo(n) = yo(ia) + dble(ii)*cell(2,1) + dble(jj)*cell(2,2) + dble(kk)*cell(2,3)
            zo(n) = zo(ia) + dble(ii)*cell(3,1) + dble(jj)*cell(3,2) + dble(kk)*cell(3,3)
            atom_name(n) = atom_name(ia)
            atom_type(n) = atom_type(ia)
!            if (luselabels) then
!              catom_label(n) = catom_label(ia)
!            else
!              atom_label(n) = atom_label(ia)
!            end if
            atom_label(n) = atom_label(ia)
            reference_number(n) = reference_number(ia)
            reference_cell(n,1) = i - 1
            reference_cell(n,2) = j - 1
            reference_cell(n,3) = k - 1
          end do
        end do loopk
      end do
    end do

    do i = 1,3
      do j = 1,3
        cell(i,j) = cell(i,j)*dble(ncell(j))
      end do
    end do
    if (ldiag) then
      write(main_output,*)
      write(main_output,'(a)') 'Output from expand_structure'
      write(main_output,'(a)') ' === Expanded structure ==='
      write(main_output,'(a)') ' New lattice vectors'
      write(main_output,*) cell(:,1)
      write(main_output,*) cell(:,2)
      write(main_output,*) cell(:,3)
    end if
    a = a*dble(ncell(1))
    b = b*dble(ncell(2))
    c = c*dble(ncell(3))
    do i = 1,newatoms
      xf(i) = xf(i)/dble(ncell(1))
      yf(i) = yf(i)/dble(ncell(2))
      zf(i) = zf(i)/dble(ncell(3))
      if (ldiag) then
!        if (luselabels) then
!         write(main_output,'(i0,2x,a4,2x,a10,2x,6f10.5,4i4)') i,atom_name(i),catom_label(i),  &
!                  xf(i),yf(i),zf(i),xo(i),yo(i),zo(i), &
!                  reference_number(i),reference_cell(i,:)
!        else
!         write(main_output,'(i0,2x,a4,2x,i0,6f10.5,4i4)') i,atom_name(i),atom_label(i),  &
!                  xf(i),yf(i),zf(i),xo(i),yo(i),zo(i), &
!                  reference_number(i),reference_cell(i,:)
!        end if
       write(main_output,'(i0,2x,a4,2x,i0,6f10.5,4i4)') i,atom_name(i),atom_label(i),  &
                  xf(i),yf(i),zf(i),xo(i),yo(i),zo(i),reference_number(i),reference_cell(i,:)
      end if
    end do

    natoms = newatoms
    if (n.ne.natoms) then
      write(6,'(a)') 'Something is wrong with the expansion of the cell'
      write(6,'(a,i0,a,i0)') 'Expected number of atoms = ',natoms,', but actual number of atoms = ',n
      stop 
    endif

    if (lsort) call sort_structure
    
    if (lreduce) ncell = ncellt

    return

  end subroutine expand_structure

!===============================================================================

  subroutine rotate_structure
!   ==============================
!--------------------------------------------------------------------------
!
!     This subroutine rotates the expanded configuration about three chosen
!     angles. The nanoparticles will be then chopped from the rotated 
!     expanded configuration
!--------------------------------------------------------------------------
    use structure_data
    use arguments
    use channel_numbers
    
    implicit none

    
    double precision :: euler(3,3), rcentre(3)
    double precision, allocatable :: xtemp(:),ytemp(:),ztemp(:)
    integer :: i
    
    
    
    euler(1,1) = cos(psi)*cos(fi)-cos(csi)*sin(fi)*sin(psi)    !These are the elements of the euler rotation matrix
    euler(1,2) = cos(psi)*sin(fi)+cos(csi)*cos(fi)*sin(psi) !Fi is the angle of rotation about the z axis
    euler(1,3) = sin(psi)*sin(csi)                            !csi about the former x axis
    euler(2,1) = -sin(psi)*cos(fi)-cos(csi)*sin(fi)*cos(psi)   !Psi about former z axis 
    euler(2,2) = -sin(psi)*sin(fi)+cos(csi)*cos(fi)*cos(psi)   !They are provided by the user in square brackets after -rotate
    euler(2,3) = cos(psi)*sin(csi)
    euler(3,1) = sin(csi)*sin(fi)
    euler(3,2) = -sin(csi)*cos(fi)
    euler(3,3) = cos(csi)
    
    
    !allocate (xfrot(natoms),yfrot(natoms),zfrot(natoms))
    allocate (xtemp(natoms),ytemp(natoms),ztemp(natoms))
    allocate (xrot(natoms),yrot(natoms),zrot(natoms))
    
    !write(6,*) cell
    
    rcentre(1) = (cell(1,1) + cell(2,1) + cell(3,1))/2.0d0
    rcentre(2) = (cell(1,2) + cell(2,2) + cell(3,2))/2.0d0
    rcentre(3) = (cell(1,3) + cell(2,3) + cell(3,3))/2.0d0

    xtemp = xo - rcentre(1)
    ytemp = yo - rcentre(2)
    ztemp = zo - rcentre(3)
    
    !translated_coordinates_loop: do i=1,natoms
    ! xtemp(i)=xo(i)-(euler(1,1)*rcentre(1)+euler(1,2)*rcentre(2)+euler(1,3)*rcentre(3))!(cell(1,1) + cell(2,1) + cell(3,1))/2.0d0!0.5*ncell(1)*cell(1,1)
     !ytemp(i)=yo(i)-(euler(2,1)*rcentre(1)+euler(2,2)*rcentre(2)+euler(2,3)*rcentre(3))!(cell(1,2) + cell(2,2) + cell(3,2))/2.0d0!0.5*ncell(2)*sqrt(cell(2,2)**2+cell(3,2)**2)
     !ztemp(i)=zo(i)-(euler(3,1)*rcentre(1)+euler(3,2)*rcentre(2)+euler(3,3)*rcentre(3))!(cell(1,3) + cell(2,3) + cell(3,3))/2.0d0!0.5*ncell(3)*cell(3,3)
    !end do translated_coordinates_loop
    
    !rotate_fractional_loop: do i=1,natoms
     !xfrot(i)=euler(1,1)*xf(i)+euler(1,2)*yf(i)+euler(1,3)*zf(i)
     !yfrot(i)=euler(2,1)*xf(i)+euler(2,2)*yf(i)+euler(2,3)*zf(i)
     !zfrot(i)=euler(3,1)*xf(i)+euler(3,2)*yf(i)+euler(3,3)*zf(i)
    !end do rotate_fractional_loop

    xrot = euler(1,1)*xtemp + euler(1,2)*ytemp + euler(1,3)*ztemp + rcentre(1)
    yrot = euler(2,1)*xtemp + euler(2,2)*ytemp + euler(2,3)*ztemp + rcentre(2)
    zrot = euler(3,1)*xtemp + euler(3,2)*ytemp + euler(3,3)*ztemp + rcentre(3)
    
    !write(6,*) sum(xo)/natoms,sum(yo)/natoms,sum(zo)/natoms
    !write(6,*) sum(xrot)/natoms,sum(yrot)/natoms,sum(zrot)/natoms
    !write(6,*)
    !write(6,*) euler(1,:)
    !write(6,*) euler(2,:)
    !write(6,*) euler(3,:)
    !write(6,*) euler(1,1)*euler(2,2)*euler(3,3) + euler(2,1)*euler(3,2)*euler(1,3) &
            ! + euler(3,1)*euler(1,2)*euler(2,3) - euler(1,1)*euler(3,2)*euler(2,3) &
            ! - euler(2,2)*euler(1,3)*euler(3,1) - euler(3,3)*euler(1,2)*euler(2,1)
    
    
    !rotate_loop: do i=1,natoms
     !xrot(i)=euler(1,1)*xo(i)+euler(1,2)*yo(i)+euler(1,3)*zo(i)+rcentre(1)-&
     !euler(1,1)*rcentre(1)-euler(1,2)*rcentre(2)-euler(1,3)*rcentre(3)!0.5*ncell(1)*cell(1,1)
     
     !yrot(i)=euler(2,1)*xo(i)+euler(2,2)*yo(i)+euler(2,3)*zo(i)+rcentre(2)-&
     !euler(2,1)*rcentre(1)-euler(2,2)*rcentre(2)-euler(2,3)*rcentre(3)!0.5*ncell(2)*sqrt(cell(2,2)**2+cell(3,2)**2)
     
     !zrot(i)=euler(3,1)*xo(i)+euler(3,2)*yo(i)+euler(3,3)*zo(i)+rcentre(3)-&
     !euler(3,1)*rcentre(1)-euler(3,2)*rcentre(2)-euler(3,3)*rcentre(3)!0.5*ncell(3)*cell(3,3)
    !end do rotate_loop
    
    !write(6,'(a,i0)') "size of xrot is",size(xrot)
    
    if (ldiag) then
      write(main_output,'(/a)') 'Output from rotate_structure'
      write(main_output,'(a/)') '===================================='
      write(main_output,'(/a/)') 'Coordinates xrot yrot zrot'
!      do i = 1,natoms
!        write(main_output,'(a)') trim(caxyz(i))
!      end do
      write(main_output,'(/a/)') 'Element, label, rotated coordinates'
      do i = 1,natoms
        write(main_output,'(a,3f12.4)') element(atom_type(i)),xrot(i),yrot(i),zrot(i)
      end do
    end if
    
    return
      
 end subroutine rotate_structure

!===============================================================================

    subroutine orient_molecules
!   ===========================

!--------------------------------------------------------------------------
!
!     This subroutine orients molecules about a random axis by different shapes, including    
!     orthorhgonal, random and flip.
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use channel_numbers
    use utilities
    
    implicit none

    integer :: ierror, i, n, m, nrigid_bodies, j
    character(200) :: text,buffer
    double precision, allocatable :: xtemp(:),ytemp(:),ztemp(:),xrotn(:),yrotn(:),zrotn(:),&
                                     xr(:), yr(:), zr(:), xcf(:), ycf(:), zcf(:)
    double precision :: r1, xq(3), rotate(3,3),volume,rvolume,rcell(3,3),u,weight,dx,dy,dz,dxo,dyo,dzo
    
    integer, allocatable :: natoms_in_rigid_body(:),rigid_body_atoms(:),atom_in_rigid_body(:,:)
    character(len=4000),allocatable :: crigid(:)
    logical :: exist
    double precision :: qa,qb,qc,qd,dorth
    
!   Check if the molecule file exists and if not, call locate_molecules
    inquire(file='molecules.dat', exist = exist)
    if (exist) then
      continue
    else
      call locate_molecules
    end if 
    
    
    open(imolecule, file=trim(wfolder)//'molecules.dat')
    ierror = 0 

    
!   Read in the number of molecules   
    read(imolecule, '(a)', iostat = ierror) buffer
    if (index(buffer, 'Number of molecules =') /= 0) then
      n = index(buffer, '=') + 1
      read(buffer(n:), *) nrigid_bodies
    end if   
    
     if (allocated(natoms_in_rigid_body)) deallocate(natoms_in_rigid_body)
    allocate(natoms_in_rigid_body(nrigid_bodies))
    allocate(xr(nrigid_bodies), yr(nrigid_bodies), zr(nrigid_bodies),crigid(nrigid_bodies)) 
    allocate(xcf(nrigid_bodies), ycf(nrigid_bodies), zcf(nrigid_bodies))
    if (allocated(xco).and.allocated(yco).and.allocated(zco)) then
           continue
        else
           allocate (xco(nrigid_bodies),yco(nrigid_bodies),zco(nrigid_bodies))
    end if
   
!   Read in the atom number within each molecule and the coordination of centre of mass
    ierror = 0
    i = 0
    do while (ierror == 0) 
      read(imolecule,'(a)',iostat=ierror) buffer
      buffer = adjustl(buffer)
      if(ierror/=0) exit
      if(len_trim(buffer)==0) cycle
       
      if ((ichar(buffer(1:1))>=ichar('0')).and.(ichar(buffer(1:1))<=ichar('9'))) then

      i = i + 1

      read(buffer,*) natoms_in_rigid_body(i),xco(i),yco(i),zco(i)
      buffer = adjustl(buffer)
        n = index(buffer,' ') ; buffer = adjustl(buffer(n:)) !pass natoms_in_rigid_body(i) 
        n = index(buffer,' ') ; buffer = adjustl(buffer(n:)) !pass xco 
        n = index(buffer,' ') ; buffer = adjustl(buffer(n:)) !pass yco   
        n = index(buffer,' ') ; buffer = adjustl(buffer(n:)) !pass zco  
      crigid(i)=buffer !now crigid(i) contains the atom number for each atom in molecule(i)
      end if
    end do
    close(imolecule)
    
    allocate(atom_in_rigid_body(nrigid_bodies,maxval(natoms_in_rigid_body))) 
    
    do i = 1, nrigid_bodies
      read(crigid(i),*) atom_in_rigid_body(i,:) !now construct a matrix which contains all the atom number for all the molecules     
    end do   

    !   From the centre of mass recalculate the xo values of each atom for each molecule
    allocate (xtemp(natoms),ytemp(natoms),ztemp(natoms))
    allocate (xrotn(natoms),yrotn(natoms),zrotn(natoms))

     volume = cell(1,1)*cell(2,2)*cell(3,3) + cell(1,2)*cell(2,3)*cell(3,1) + &
				cell(1,3)*cell(2,1)*cell(3,2) - cell(1,1)*cell(2,3)*cell(3,2) - &
				cell(1,3)*cell(2,2)*cell(3,1) - cell(1,2)*cell(2,1)*cell(3,3)
  
    ! Calculate the Cartesian coordinates of centre of mass
	   rcell(1,1) = cell(2,2)*cell(3,3) - cell(2,3)*cell(3,2) ! correct
	   rcell(1,2) = cell(2,3)*cell(3,1) - cell(2,1)*cell(3,3) ! correct
	   rcell(1,3) = cell(2,1)*cell(3,2) - cell(2,2)*cell(3,1) ! correct
	   rcell(2,1) = cell(3,2)*cell(1,3) - cell(3,3)*cell(1,2) ! correct
	   rcell(2,2) = cell(3,3)*cell(1,1) - cell(3,1)*cell(1,3) ! correct
	   rcell(2,3) = cell(3,1)*cell(1,2) - cell(3,2)*cell(1,1) ! correct
	   rcell(3,1) = cell(1,2)*cell(2,3) - cell(1,3)*cell(2,2) ! correct
	   rcell(3,2) = cell(1,3)*cell(2,1) - cell(1,1)*cell(2,3) ! correct
	   rcell(3,3) = cell(1,1)*cell(2,2) - cell(1,2)*cell(2,1) ! correct

	   rcell = rcell/volume

     do i = 1, nrigid_bodies
       xcf(i) = xf(atom_in_rigid_body(i,1))
       ycf(i) = yf(atom_in_rigid_body(i,1))
       zcf(i) = zf(atom_in_rigid_body(i,1))         
     end do
     
     do i = 1, nrigid_bodies 
!       weight = element_mass(atom_type(atom_in_rigid_body(i,1)))
       do m=2, natoms_in_rigid_body(i)
            j = atom_in_rigid_body(i,m)
            dx = xcf(i)/(m-1) - xf(j) + 1.5d0
            dx = dx - aint(dx) - 0.5d0
            dy = ycf(i)/(m-1) - yf(j) + 1.5d0
            dy = dy - aint(dy) - 0.5d0
            dz = zcf(i)/(m-1) - zf(j) + 1.5d0
            dz = dz - aint(dz) - 0.5d0
            write(506,*) dx,dy,dz
            xf(j) = xcf(i)/(m-1) - dx
            yf(j) = ycf(i)/(m-1) - dy
            zf(j) = zcf(i)/(m-1) - dz
            xcf(i) = xcf(i) + xf(j)
            ycf(i) = ycf(i) + yf(j)
            zcf(i) = zcf(i) + zf(j)
 !          xcf(i) = xcf(i)*weight + xf(j)*element_mass(atom_type(j))
 !          ycf(i) = ycf(i)*weight + yf(j)*element_mass(atom_type(j))
 !          zcf(i) = zcf(i)*weight + zf(j)*element_mass(atom_type(j))
 !          weight = weight + element_mass(atom_type(j))
 !          xcf(i) = xcf(i)/weight 
 !          ycf(i) = ycf(i)/weight 
 !          zcf(i) = zcf(i)/weight 
       end do
        xcf(i) = xcf(i)/dble(natoms_in_rigid_body(i))
        ycf(i) = ycf(i)/dble(natoms_in_rigid_body(i))
        zcf(i) = zcf(i)/dble(natoms_in_rigid_body(i))
     end do 
     
    xco = cell(1,1)*xcf + cell(1,2)*ycf + cell(1,3)*zcf
    yco = cell(2,1)*xcf + cell(2,2)*ycf + cell(2,3)*zcf
    zco = cell(3,1)*xcf + cell(3,2)*ycf + cell(3,3)*zcf

    xo = cell(1,1)*xf + cell(1,2)*yf + cell(1,3)*zf
    yo = cell(2,1)*xf + cell(2,2)*yf + cell(2,3)*zf
    zo = cell(3,1)*xf + cell(3,2)*yf + cell(3,3)*zf

do i = 1,nrigid_bodies
  write(501,'(a2,7x,a2,i0,4x,3f12.6)') 'Xe','Xe',i,xcf(i),ycf(i),zcf(i)
end do
      
     do i = 1, nrigid_bodies
       do m=1, natoms_in_rigid_body(i) 
           xtemp(atom_in_rigid_body(i,m)) = xo(atom_in_rigid_body(i,m)) - xco(i)
           ytemp(atom_in_rigid_body(i,m)) = yo(atom_in_rigid_body(i,m)) - yco(i)
           ztemp(atom_in_rigid_body(i,m)) = zo(atom_in_rigid_body(i,m)) - zco(i)
    write(600,*)  SQRT(xtemp(atom_in_rigid_body(i,m))**2 + ytemp(atom_in_rigid_body(i,m))**2 + &
     ztemp(atom_in_rigid_body(i,m))**2)
       end do
     end do
   
       
    do i = 1, nrigid_bodies
      select case(trim(oshape))
        case('orth','90')
          call random_number(xq)
          xq = xq/maxval(xq)
          xq = int(xq)
          call random_number(r1)
          r1 = r1*4.0d0
          dorth = dasin(1.0d0)*int(r1)
      qa = cos(dorth/2); qb = sin(dorth/2)*xq(1)
      qc = sin(dorth/2)*xq(2); qd = sin(dorth/2)*xq(3)
        case('random')
          call random_number(xq)
          xq = 2.0d0*xq - 1.0d0
          u = sqrt(xq(1)**2 + xq(2)**2 + xq(3)**2)
          xq = xq/u
          call random_number(r1)
          r1 = 2.0d0*r1 - 1.0d0
          dorth = asin(r1)
          qa = random_normal()
          qb = random_normal()
          qc = random_normal() 
          qd = random_normal()
          u = sqrt(qa**2 + qb**2 + qc**2 + qd**2)
          qa = qa/u
          qb = qb/u
          qc = qc/u
          qd = qd/u
        case('flipx','180x')
          call random_number(r1)
          dorth = dble(nint(r1))*acos(-1.0d0)
          xq = (/ 1.0d0, 0.0d0, 0.0d0 /)
          qa = cos(dorth/2); qb = sin(dorth/2)*xq(1)
          qc = sin(dorth/2)*xq(2); qd = sin(dorth/2)*xq(3)
        case('flipy','180y')
          call random_number(r1)
          dorth = dble(nint(r1))*acos(-1.0d0)
          xq = (/ 0.0d0, 1.0d0, 0.0d0 /)
          qa = cos(dorth/2); qb = sin(dorth/2)*xq(1)
          qc = sin(dorth/2)*xq(2); qd = sin(dorth/2)*xq(3)
        case('flipz','180z')
          call random_number(r1)
          dorth = dble(nint(r1))*acos(-1.0d0)
          xq = (/ 0.0d0, 0.0d0, 1.0d0 /)
          qa = cos(dorth/2); qb = sin(dorth/2)*xq(1)
          qc = sin(dorth/2)*xq(2); qd = sin(dorth/2)*xq(3)
        case('flipr','180r')
          call random_number(xq)
          xq = 2.0d0*xq - 1.0d0
          u = sqrt(xq(1)**2 + xq(2)**2 + xq(3)**2)
          xq = xq/u
          call random_number(r1)
          r1 = 2.0d0*r1 - 1.0d0
          dorth = acos(-1.0d0)
          qa = cos(dorth/2); qb = sin(dorth/2)*xq(1)
          qc = sin(dorth/2)*xq(2); qd = sin(dorth/2)*xq(3)
      end select 

      rotate(1,1) = qa**2+qb**2-qc**2-qd**2  ; rotate(1,2) = 2*qb*qc-2*qa*qd  ; rotate(1,3) = 2*qb*qd+2*qa*qc;       
      rotate(2,1) = 2*qb*qc+2*qa*qd  ;  rotate(2,2) = qa**2-qb**2+qc**2-qd**2  ; rotate(2,3) = 2*qc*qd-2*qa*qb
      rotate(3,1) = 2*qb*qd-2*qa*qc  ; rotate(3,2) = 2*qc*qd+2*qa*qb  ; rotate(3,3) = qa**2-qb**2-qc**2+qd**2
      
!         do m = 1, natoms_in_rigid_body(i) 
!           xrotn(atom_in_rigid_body(i,m)) = xtemp(atom_in_rigid_body(i,m))*rotate(1,1)+&
!                   ytemp(atom_in_rigid_body(i,m))*rotate(1,2)+&
!                   ztemp(atom_in_rigid_body(i,m))*rotate(1,3)+xco(i)
!           yrotn(atom_in_rigid_body(i,m)) = xtemp(atom_in_rigid_body(i,m))*rotate(2,1)+&
!                   ytemp(atom_in_rigid_body(i,m))*rotate(2,2)+&
!                   ztemp(atom_in_rigid_body(i,m))*rotate(2,3)+yco(i)
!           zrotn(atom_in_rigid_body(i,m)) = xtemp(atom_in_rigid_body(i,m))*rotate(3,1)+&
!                   ytemp(atom_in_rigid_body(i,m))*rotate(3,2)+&
!                   ztemp(atom_in_rigid_body(i,m))*rotate(3,3)+zco(i)
!         end do
        do m = 1, natoms_in_rigid_body(i) 
          xo(atom_in_rigid_body(i,m)) = xtemp(atom_in_rigid_body(i,m))*rotate(1,1)+&
                  ytemp(atom_in_rigid_body(i,m))*rotate(1,2)+&
                  ztemp(atom_in_rigid_body(i,m))*rotate(1,3)+xco(i)
          yo(atom_in_rigid_body(i,m)) = xtemp(atom_in_rigid_body(i,m))*rotate(2,1)+&
                  ytemp(atom_in_rigid_body(i,m))*rotate(2,2)+&
                  ztemp(atom_in_rigid_body(i,m))*rotate(2,3)+yco(i)
          zo(atom_in_rigid_body(i,m)) = xtemp(atom_in_rigid_body(i,m))*rotate(3,1)+&
                  ytemp(atom_in_rigid_body(i,m))*rotate(3,2)+&
                  ztemp(atom_in_rigid_body(i,m))*rotate(3,3)+zco(i)
        end do
    end do
    
      volume = cell(1,1)*cell(2,2)*cell(3,3) + cell(1,2)*cell(2,3)*cell(3,1) + &
				cell(1,3)*cell(2,1)*cell(3,2) - cell(1,1)*cell(2,3)*cell(3,2) - &
				cell(1,3)*cell(2,2)*cell(3,1) - cell(1,2)*cell(2,1)*cell(3,3)

	   rcell(1,1) = cell(2,2)*cell(3,3) - cell(2,3)*cell(3,2) ! correct
	   rcell(1,2) = cell(2,3)*cell(3,1) - cell(2,1)*cell(3,3) ! correct
	   rcell(1,3) = cell(2,1)*cell(3,2) - cell(2,2)*cell(3,1) ! correct
	   rcell(2,1) = cell(3,2)*cell(1,3) - cell(3,3)*cell(1,2) ! correct
	   rcell(2,2) = cell(3,3)*cell(1,1) - cell(3,1)*cell(1,3) ! correct
	   rcell(2,3) = cell(3,1)*cell(1,2) - cell(3,2)*cell(1,1) ! correct
	   rcell(3,1) = cell(1,2)*cell(2,3) - cell(1,3)*cell(2,2) ! correct
	   rcell(3,2) = cell(1,3)*cell(2,1) - cell(1,1)*cell(2,3) ! correct
	   rcell(3,3) = cell(1,1)*cell(2,2) - cell(1,2)*cell(2,1) ! correct

	   rcell = rcell/volume

	   rvolume = rcell(1,1)*rcell(2,2)*rcell(3,3) + rcell(1,2)*rcell(2,3)*rcell(3,1) + &
				 rcell(1,3)*rcell(2,1)*rcell(3,2) - rcell(1,1)*rcell(2,3)*rcell(3,2) - &
				 rcell(1,3)*rcell(2,2)*rcell(3,1) - rcell(1,2)*rcell(2,1)*rcell(3,3)
    do i = 1, nrigid_bodies
       do m=1, natoms_in_rigid_body(i) 
         j = atom_in_rigid_body(i,m)
        xf(j) = rcell(1,1)*xo(j) + rcell(1,2)*yo(j) + rcell(1,3)*zo(j)
        yf(j) = rcell(2,1)*xo(j) + rcell(2,2)*yo(j) + rcell(2,3)*zo(j)
        zf(j) = rcell(3,1)*xo(j) + rcell(3,2)*yo(j) + rcell(3,3)*zo(j)
     end do
  end do
  
    
 end subroutine orient_molecules    

!===============================================================================
     
    subroutine transform_structure
!   ==============================

!--------------------------------------------------------------------------
!
!     This subroutine transforms the basic structure to a new cell of cell
!     vectors. Unlike the supercell approach, this allows for rotation of
!     axes and the like; in the end it should supercede expand_structure.
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use channel_numbers

    implicit none

    integer :: i,j,k,ii,jj,kk,n,ia,newatoms,mvolume,nmax
    integer, allocatable :: atom_typet(:),atom_labelt(:),rnumt(:),rcellt(:,:)
    double precision,allocatable :: xft(:),yft(:),zft(:),xot(:),yot(:),zot(:)
    character(len=10), allocatable :: catom_labelt(:)
    character(len=4), allocatable :: atom_namet(:)
    double precision :: rtransform(3,3),xyzp(3)

    allocate(xft(natoms),yft(natoms),zft(natoms))
    allocate(xot(natoms),yot(natoms),zot(natoms))
    allocate(atom_namet(natoms))
    allocate(atom_typet(natoms))
!    if (luselabels) then
!      allocate(catom_labelt(natoms))
!    else
!      allocate(atom_labelt(natoms))
!    end if
    allocate(atom_labelt(natoms))
    allocate(rnumt(natoms),rcellt(natoms,3))
    xft = xf(1:natoms)
    yft = yf(1:natoms)
    zft = zf(1:natoms)
    xot = xo(1:natoms)
    yot = yo(1:natoms)
    zot = zo(1:natoms)
    atom_namet = atom_name(1:natoms)
    atom_typet = atom_type(1:natoms)
!    if (luselabels) then
!      catom_labelt = catom_label(1:natoms)
!    else
!      atom_labelt = atom_label(1:natoms)
!    end if
    atom_labelt = atom_label(1:natoms)
    rnumt = reference_number(1:natoms)
    rcellt = reference_cell(1:natoms,:)
    if (allocated(xf)) deallocate(xf,yf,zf)
    if (allocated(xo)) deallocate(xo,yo,zo)
    if (allocated(atom_name)) deallocate(atom_name)
    if (allocated(atom_type)) deallocate(atom_type)
    if (allocated(atom_label)) deallocate(atom_label)
    if (allocated(catom_label)) deallocate(catom_label)
    if (allocated(reference_number)) deallocate(reference_number)
    if (allocated(reference_cell)) deallocate(reference_cell)

    mvolume = mtransform(1,1)*mtransform(2,2)*mtransform(3,3) + &
              mtransform(1,2)*mtransform(2,3)*mtransform(3,1) + &
              mtransform(2,1)*mtransform(1,3)*mtransform(3,2) - &
              mtransform(1,1)*mtransform(2,3)*mtransform(3,2) - &
              mtransform(2,2)*mtransform(1,3)*mtransform(3,1) - &
              mtransform(3,3)*mtransform(1,2)*mtransform(2,1)
    newatoms = natoms*mvolume
    rtransform(1,1) = dble(mtransform(2,2)*mtransform(3,3)-mtransform(2,3)*mtransform(3,2))
    rtransform(2,2) = dble(mtransform(1,1)*mtransform(3,3)-mtransform(1,3)*mtransform(3,1))
    rtransform(3,3) = dble(mtransform(1,1)*mtransform(2,2)-mtransform(1,2)*mtransform(2,1))
    rtransform(2,1) = -dble(mtransform(2,1)*mtransform(3,3)-mtransform(2,3)*mtransform(3,1))
    rtransform(3,1) = dble(mtransform(2,1)*mtransform(3,2)-mtransform(2,2)*mtransform(3,1))
    rtransform(1,2) = -dble(mtransform(1,2)*mtransform(3,3)-mtransform(1,3)*mtransform(3,2))
    rtransform(3,2) = -dble(mtransform(1,1)*mtransform(3,2)-mtransform(1,2)*mtransform(3,1))
    rtransform(1,3) = dble(mtransform(1,2)*mtransform(2,3)-mtransform(1,3)*mtransform(2,2))
    rtransform(2,3) = -dble(mtransform(1,1)*mtransform(2,3)-mtransform(1,3)*mtransform(2,1))
    rtransform = rtransform/dble(mvolume)
    write(6,*) mvolume
    write(6,'(3i7)') mtransform
    write(6,*)
    write(6,'(3f7.2)') rtransform
    write(6,*)
    write(6,'(3f7.2)') matmul(mtransform,rtransform)
    write(6,*)
    write(6,'(3f7.2)') cell
    cell = matmul(mtransform,cell)
    write(6,*)
    write(6,'(3f7.2)') cell
    a = dsqrt(cell(1,1)**2 + cell(2,1)**2 + cell(3,1)**2)
    b = dsqrt(cell(1,2)**2 + cell(2,2)**2 + cell(3,2)**2)
    c = dsqrt(cell(1,3)**2 + cell(2,3)**2 + cell(3,3)**2)
    alpha = acos((cell(1,2)*cell(1,3) + cell(2,2)*cell(2,3) + cell(3,2)*cell(3,3))/(b*c))
    beta  = acos((cell(1,1)*cell(1,3) + cell(2,1)*cell(2,3) + cell(3,1)*cell(3,3))/(a*c))
    gamma = acos((cell(1,1)*cell(1,2) + cell(2,1)*cell(2,2) + cell(3,1)*cell(3,2))/(a*b))
    alpha = alpha*180.0d0/acos(-1.0d0)
    beta = beta*180.0d0/acos(-1.0d0)
    gamma = gamma*180.0d0/acos(-1.0d0)
    write(6,'(3f8.3)') a,b,c
    write(6,'(3f8.3)') alpha,beta,gamma
    nmax = max(maxval(mtransform),abs(minval(mtransform)))
    write(6,*) nmax
    do i = -nmax,nmax
      do j = -nmax,nmax
        do k = -nmax,nmax
        end do
      end do
    end do

stop

    do i = -1,nelements
      n_elements(i) = n_elements(i)*ncell(1)*ncell(2)*ncell(3)
    end do
    allocate(xf(newatoms),yf(newatoms),zf(newatoms))
    allocate(xo(newatoms),yo(newatoms),zo(newatoms))
    allocate(atom_name(newatoms))
    allocate(atom_type(newatoms))
!    if (luselabels) then
!      allocate(catom_label(newatoms))
!    else
!      allocate(atom_label(newatoms))
!    end if
    allocate(atom_label(newatoms))
    allocate(reference_number(newatoms))
    allocate(reference_cell(newatoms,3))
    xf(1:natoms) = xft
    yf(1:natoms) = yft
    zf(1:natoms) = zft
    xo(1:natoms) = xot
    yo(1:natoms) = yot
    zo(1:natoms) = zot
    atom_name(1:natoms) = atom_namet
    atom_type(1:natoms) = atom_typet
!    if (luselabels) then
!      catom_label(1:natoms) = catom_labelt
!    else
!      atom_label(1:natoms) = atom_labelt
!    end if
    atom_label(1:natoms) = atom_labelt
    reference_number(1:natoms) = rnumt
    reference_cell(1:natoms,:) = rcellt
    if (allocated(xft)) deallocate(xft,yft,zft)
    if (allocated(xot)) deallocate(xot,yot,zot)
    if (allocated(atom_namet)) deallocate(atom_namet)
    if (allocated(atom_typet)) deallocate(atom_typet)
    if (allocated(atom_labelt)) deallocate(atom_labelt)
    if (allocated(catom_labelt)) deallocate(catom_labelt)
    if (allocated(rnumt)) deallocate(rnumt)
    if (allocated(rcellt)) deallocate(rcellt)

   do ia = 1,natoms
     xyzp(1) = xf(ia) ; xyzp(2) = yf(ia) ; xyzp(3) = zf(ia)
     write(6,'(3f8.4)') xyzp
     xyzp = matmul(transpose(rtransform),xyzp)
     write(6,'(3f8.4)') xyzp
     write(6,*)
  end do
  stop


    n = natoms
    do i = 1,ncell(1)
    ii = i-1
      do j = 1,ncell(2)
      jj = j-1
        loopk: do k = 1,ncell(3)
          if ((i.eq.1).and.(j.eq.1).and.(k.eq.1)) cycle loopk
          kk = k-1
          do ia = 1,natoms
            n = n + 1
            xf(n) = xf(ia) + dble(ii)
            yf(n) = yf(ia) + dble(jj)
            zf(n) = zf(ia) + dble(kk)
            xo(n) = xo(ia) + dble(ii)*cell(1,1) + dble(jj)*cell(1,2) + dble(kk)*cell(1,3)
            yo(n) = yo(ia) + dble(ii)*cell(2,1) + dble(jj)*cell(2,2) + dble(kk)*cell(2,3)
            zo(n) = zo(ia) + dble(ii)*cell(3,1) + dble(jj)*cell(3,2) + dble(kk)*cell(3,3)
            atom_name(n) = atom_name(ia)
            atom_type(n) = atom_type(ia)
!            if (luselabels) then
!              catom_label(n) = catom_label(ia)
!            else
!              atom_label(n) = atom_label(ia)
!            end if
            atom_label(n) = atom_label(ia)
            reference_number(n) = reference_number(ia)
            reference_cell(n,1) = i - 1
            reference_cell(n,2) = j - 1
            reference_cell(n,3) = k - 1
          end do
        end do loopk
      end do
    end do

    do i = 1,3
      do j = 1,3
        cell(i,j) = cell(i,j)*dble(ncell(j))
      end do
    end do
    if (ldiag) then
      write(main_output,*)
      write(main_output,'(a)') 'Output from expand_structure'
      write(main_output,'(a)') ' === Expanded structure ==='
      write(main_output,'(a)') ' New lattice vectors'
      write(main_output,*) cell(:,1)
      write(main_output,*) cell(:,2)
      write(main_output,*) cell(:,3)
    end if
    a = a*dble(ncell(1))
    b = b*dble(ncell(2))
    c = c*dble(ncell(3))
    do i = 1,newatoms
      xf(i) = xf(i)/dble(ncell(1))
      yf(i) = yf(i)/dble(ncell(2))
      zf(i) = zf(i)/dble(ncell(3))
      if (ldiag) &
         write(main_output,'(i5,2x,a4,2x,6f10.5,4i4)') i,atom_name(i),  &
                  xf(i),yf(i),zf(i),xo(i),yo(i),zo(i), &
                  reference_number(i),reference_cell(i,:)
    end do

    natoms = newatoms
    if (n.ne.natoms) then
      if (ldiag) write(main_output,*) ' **** Something is wrong with the expansion of the cell'
      stop 'Something is wrong with the expansion of the cell'
    endif

    if (lsort) call sort_structure

    return

  end subroutine transform_structure

!===============================================================================

    subroutine allocate_grid
!   ========================

!--------------------------------------------------------------------------
!
!     This subroutine divides the structure into grid cells and allocates
!     each atom to a cell, then it creates a neighbour list for each atom
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use channel_numbers

    implicit none

    integer :: i,j,jg,n,iatom,nx,ny,nz,ix,iy,iz,jx,jy,jz,kx,ky,kz,nmax
    integer, allocatable :: ngatoms(:,:,:),neighsmax(:,:,:)
    double precision :: dx,dxo,dy,dyo,dz,dzo,r

    if (ldiag) then
      write(main_output,'(/a)') 'Output from allocate grid'
      write(main_output,'(a/)') '========================='
      write(main_output,'(a,f7.3)') 'Grid cell size = ',gridcellsize
    end if

   if(allocated(gridcell)) deallocate(gridcell)
   if(allocated(ngatoms)) deallocate(ngatoms)
   if(allocated(neighsmax)) deallocate(neighsmax)
   if(allocated(gridatoms)) deallocate(gridatoms)
   if(allocated(neighbours)) deallocate(neighbours)

   allocate(gridcell(natoms,3))

   if (ladjust_xyz_signs) then
     do i = 1,natoms
       do while (xf(i)<0.0d0)
         xf(i) = xf(i) + 1.0d0
       end do
       do while (yf(i)<0.0d0)
         yf(i) = yf(i) + 1.0d0
       end do
       do while (zf(i)<0.0d0)
         zf(i) = zf(i) + 1.0d0
       end do
       do while (xf(i)>=1.0d0)
         xf(i) = xf(i) - 1.0d0
       end do
       do while (yf(i)>=1.0d0)
         yf(i) = yf(i) - 1.0d0
       end do
       do while (zf(i)>=1.0d0)
         zf(i) = zf(i) - 1.0d0
       end do
     end do
   end if

    nx = max(int(a/gridcellsize),3)
    ny = max(int(b/gridcellsize),3)
    nz = max(int(c/gridcellsize),3)
    allocate(ngatoms(nx,ny,nz))
!    write(6,'(a,i0,a,i0,a,i0)') 'Grid shape = ',nx,'x',ny,'x',nz

!   Allocate grid cell against each atom
    do iatom = 1,natoms
      gridcell(iatom,1) = int(xf(iatom)*dble(nx)) + 1
      gridcell(iatom,2) = int(yf(iatom)*dble(ny)) + 1
      gridcell(iatom,3) = int(zf(iatom)*dble(nz)) + 1
      if (gridcell(iatom,1)<0) gridcell(iatom,1) = gridcell(iatom,1) + nx
      if (gridcell(iatom,2)<0) gridcell(iatom,2) = gridcell(iatom,2) + ny
      if (gridcell(iatom,3)<0) gridcell(iatom,3) = gridcell(iatom,3) + nz
!write(6,'(i4,3f6.3,3i4)') iatom,xf(iatom),yf(iatom),zf(iatom),gridcell(iatom,:)
    end do

!   count number of atoms within each grid cell
    ngatoms = 0
    do iatom = 1,natoms
      ix = gridcell(iatom,1) ; iy = gridcell(iatom,2) ; iz = gridcell(iatom,3)
      ngatoms(ix,iy,iz) = ngatoms(ix,iy,iz) + 1
    end do

!   Assign atoms to each cell
    nmax = maxval(ngatoms)
    allocate(gridatoms(nx,ny,nz,0:nmax))
    gridatoms = 0
    do iatom = 1,natoms
      ix = gridcell(iatom,1) ; iy = gridcell(iatom,2) ; iz = gridcell(iatom,3)
      gridatoms(ix,iy,iz,0) = gridatoms(ix,iy,iz,0) + 1
      i = gridatoms(ix,iy,iz,0)
      gridatoms(ix,iy,iz,i) = iatom
    end do

!do ix = 1,nx
!do iy = 1,ny
!do iz = 1,nz
!write(6,'(a)') '-----------'
!write(6,*) ix,iy,iz
!do i = 1, gridatoms(ix,iy,iz,0)
!iatom = gridatoms(ix,iy,iz,i)
!write(6,'(i3,3f7.3)') iatom,xf(iatom),yf(iatom),zf(iatom)
!end do
!end do
!end do
!end do

allocate(neighsmax(nx,ny,nz))
neighsmax = 0
do ix = 1,nx
do iy = 1,ny
do iz = 1,nz
     do jx = -1,1
       kx = modulo(ix+jx-1,nx) + 1
       do jy = -1,1
         ky = modulo(iy+jy-1,ny) + 1
         do jz = -1,1
           kz = modulo(iz+jz-1,nz) + 1
           neighsmax(ix,iy,iz) = neighsmax(ix,iy,iz) + gridatoms(kx,ky,kz,0)
         end do
       end do
     end do
end do
end do
end do

!  Count neighbours for each atom
   nmax = 0
   do i = 1,natoms
     ix = gridcell(i,1) ; iy = gridcell(i,2) ; iz = gridcell(i,3)
     n = 0
     do jx = -1,1
       kx = modulo(ix+jx-1,nx) + 1
       do jy = -1,1
         ky = modulo(iy+jy-1,ny) + 1
         do jz = -1,1
           kz = modulo(iz+jz-1,nz) + 1
           do jg = 1,gridatoms(kx,ky,kz,0)
             j = gridatoms(kx,ky,kz,jg)
             if (i==j) cycle
             dx = xf(i) - xf(j) + 1.5d0
             dx = dx - aint(dx) - 0.5d0
             dy = yf(i) - yf(j) + 1.5d0
             dy = dy - aint(dy) - 0.5d0
             dz = zf(i) - zf(j) + 1.5d0
             dz = dz - aint(dz) - 0.5d0
             dxo = cell(1,1)*dx + cell(1,2)*dy + cell(1,3)*dz
             dyo = cell(2,1)*dx + cell(2,2)*dy + cell(2,3)*dz
             dzo = cell(3,1)*dx + cell(3,2)*dy + cell(3,3)*dz
             r = dsqrt(dxo**2 + dyo**2 + dzo**2)
             if (r>gridcellsize) cycle
             n = n + 1
             if (n>nmax) nmax = n
           end do
         end do
       end do
     end do
   end do
   allocate(neighbours(natoms,0:nmax))
   neighbours = 0

!  Create neighbour list for each atom
   do i = 1,natoms
     n = 0
     ix = gridcell(i,1) ; iy = gridcell(i,2) ; iz = gridcell(i,3)
     do jx = -1,1
       kx = modulo(ix+jx-1,nx) + 1
       do jy = -1,1
         ky = modulo(iy+jy-1,ny) + 1
         do jz = -1,1
           kz = modulo(iz+jz-1,nz) + 1
           do jg = 1,gridatoms(kx,ky,kz,0)
             j = gridatoms(kx,ky,kz,jg)
             if (i==j) cycle
             dx = xf(i) - xf(j) + 1.5d0
             dx = dx - aint(dx) - 0.5d0
             dy = yf(i) - yf(j) + 1.5d0
             dy = dy - aint(dy) - 0.5d0
             dz = zf(i) - zf(j) + 1.5d0
             dz = dz - aint(dz) - 0.5d0
             dxo = cell(1,1)*dx + cell(1,2)*dy + cell(1,3)*dz
             dyo = cell(2,1)*dx + cell(2,2)*dy + cell(2,3)*dz
             dzo = cell(3,1)*dx + cell(3,2)*dy + cell(3,3)*dz
             r = dsqrt(dxo**2 + dyo**2 + dzo**2)
             if (r>gridcellsize) cycle
             n = n + 1
             neighbours(i,n) = j
             neighbours(i,0) = n
           end do
         end do
       end do
     end do
   end do

!  do i = 1,natoms
!    write(6,*) element(atom_type(i)),(element(atom_type(neighbours(i,j))),j=1,neighbours(i,0))
!    do j = 1,neighbours(i,0)
!      dxo = xo(neighbours(i,j))-xo(i)
!      dyo = yo(neighbours(i,j))-yo(i)
!      dzo = zo(neighbours(i,j))-zo(i)
!      r = dsqrt(dxo**2 + dyo**2 + dzo**2)
!      write(6,*) element(atom_type(i)),element(atom_type(neighbours(i,j))),r
!    end do
!  end do

  return

  end subroutine allocate_grid

!===============================================================================

    subroutine rectangulate
!   =======================

!--------------------------------------------------------------------------
!
!     This subroutine expands the structure by converting a lattice on the
!     hexagonal setting to a C-cented orthogonal cell
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use channel_numbers

    implicit none

    integer :: i,newatoms
    integer, allocatable :: atom_typet(:),atom_labelt(:),rnumt(:),rcellt(:,:)
    double precision,allocatable :: xft(:),yft(:),zft(:),xot(:),yot(:),zot(:)
    character(len=10), allocatable :: catom_labelt(:)
    character(len=4), allocatable :: atom_namet(:)
    double precision :: cellt(3,3),transfm(3,3),pi,volume

    pi = acos(-1.0d0)

    allocate(xft(natoms),yft(natoms),zft(natoms))
    allocate(xot(natoms),yot(natoms),zot(natoms))
    allocate(atom_namet(natoms))
    allocate(atom_typet(natoms))
!    if (luselabels) then
!      allocate(catom_labelt(natoms))
!    else
!      allocate(atom_labelt(natoms))
!    end if
    allocate(atom_labelt(natoms))
    allocate(rnumt(natoms),rcellt(natoms,3))
    xft = xf(1:natoms)
    yft = yf(1:natoms)
    zft = zf(1:natoms)
    xot = xo(1:natoms)
    yot = yo(1:natoms)
    zot = zo(1:natoms)
    atom_namet = atom_name(1:natoms)
    atom_typet = atom_type(1:natoms)
!    if (luselabels) then
!      catom_labelt = catom_label(1:natoms)
!    else
!      atom_labelt = atom_label(1:natoms)
!    end if
    atom_labelt = atom_label(1:natoms)
    rnumt = reference_number(1:natoms)
    rcellt = reference_cell(1:natoms,:)
    if (allocated(xf)) deallocate(xf,yf,zf)
    if (allocated(xo)) deallocate(xo,yo,zo)
    if (allocated(atom_name)) deallocate(atom_name)
    if (allocated(atom_type)) deallocate(atom_type)
    if (allocated(atom_label)) deallocate(atom_label)
    if (allocated(catom_label)) deallocate(catom_label)
    if (allocated(reference_number)) deallocate(reference_number)
    if (allocated(reference_cell)) deallocate(reference_cell)

    newatoms = natoms*2

    do i = -1,nelements
      n_elements(i) = n_elements(i)*2
    end do

    allocate(xf(newatoms),yf(newatoms),zf(newatoms))
    allocate(xo(newatoms),yo(newatoms),zo(newatoms))
    allocate(atom_name(newatoms))
    allocate(atom_type(newatoms))
!    if (luselabels) then
!      allocate(catom_label(newatoms))
!    else
!      allocate(atom_label(newatoms))
!    end if
    allocate(atom_label(newatoms))
    allocate(reference_number(newatoms))
    allocate(reference_cell(newatoms,3))

    transfm = 0.0d0
    transfm(1,1) = 2.0d0
    transfm(1,2) = 1.0d0
    transfm(2,2) = 1.0d0
    transfm(3,3) = 1.0d0

    cellt = matmul(transfm,transpose(cell))
    cell = transpose(cellt)

    a = dsqrt(cell(1,1)**2 + cell(2,1)**2 + cell(3,1)**2)
    b = dsqrt(cell(1,2)**2 + cell(2,2)**2 + cell(3,2)**2)
    c = dsqrt(cell(1,3)**2 + cell(2,3)**2 + cell(3,3)**2)
    alpha = (180.0d0/pi)*dacos((cell(1,2)*cell(1,3) + cell(2,2)*cell(2,3) + &
                cell(3,2)*cell(3,3))/(b*c))
    beta  = (180.0d0/pi)*dacos((cell(1,1)*cell(1,3) + cell(2,1)*cell(2,3) + &
                cell(3,1)*cell(3,3))/(a*c))
    gamma = (180.0d0/pi)*dacos((cell(1,1)*cell(1,2) + cell(2,1)*cell(2,2) + &
                cell(3,1)*cell(3,2))/(a*b))
    volume = cell(1,1)*cell(2,2)*cell(3,3) + cell(2,1)*cell(3,2)*cell(1,3) + &
           cell(3,1)*cell(1,2)*cell(2,3) - cell(1,1)*cell(3,2)*cell(2,3) - &
           cell(2,1)*cell(1,2)*cell(3,3) - cell(3,1)*cell(2,2)*cell(1,3)

    transfm = 0.0d0
    transfm(1,1) = 0.5d0
    transfm(2,1) = -0.5d0
    transfm(2,2) = 1.0d0
    transfm(3,3) = 1.0d0

    do i = 1,natoms
      xf(i) = transfm(1,1)*xft(i)
      xf(i+natoms) = xf(i) + 0.5d0
      yf(i) = transfm(2,1)*xft(i) + transfm(2,2)*yft(i)
      yf(i+natoms) = yf(i) + 0.5d0
      zf(i) = zft(i)
      zf(i+natoms) = zf(i)
    end do

    do i = 1,newatoms
      if (xf(i)<0.0d0) xf(i) = xf(i) + 1.0d0
      if (xf(i)>1.0d0) xf(i) = xf(i) - 1.0d0
      if (yf(i)<0.0d0) yf(i) = yf(i) + 1.0d0
      if (yf(i)>1.0d0) yf(i) = yf(i) - 1.0d0
      if (zf(i)<0.0d0) zf(i) = zf(i) + 1.0d0
      if (zf(i)>1.0d0) zf(i) = zf(i) - 1.0d0
    end do

    do i = 1,newatoms
       xo(i) = cell(1,1)*xf(i) + cell(1,2)*yf(i) + cell(1,3)*zf(i)
       yo(i) = cell(2,1)*xf(i) + cell(2,2)*yf(i) + cell(2,3)*zf(i)
       zo(i) = cell(3,1)*xf(i) + cell(3,2)*yf(i) + cell(3,3)*zf(i)
    end do

    atom_name(1:natoms) = atom_namet
    atom_name(natoms+1:newatoms) = atom_namet
    atom_type(1:natoms) = atom_typet
    atom_type(natoms+1:newatoms) = atom_typet
!    if (luselabels) then
!      catom_label(1:natoms) = catom_labelt
!      catom_label(natoms+1:newatoms) = catom_labelt
!    else
!      atom_label(1:natoms) = atom_labelt
!      atom_label(natoms+1:newatoms) = atom_labelt
!    end if
    atom_label(1:natoms) = atom_labelt
    atom_label(natoms+1:newatoms) = atom_labelt
    reference_number(1:natoms) = rnumt
    reference_number(natoms+1:newatoms) = rnumt
    reference_cell(1:natoms,:) = rcellt
    reference_cell(natoms+1:newatoms,:) = rcellt

    natoms = newatoms

    if (ldiag) then
      write(main_output,*)
      write(main_output,'(a)') 'Output from rectangulate'
      write(main_output,'(a)') ' === Expanded structure ==='
      write(main_output,'(a)') ' New lattice vectors'
      write(main_output,*) cell(:,1)
      write(main_output,*) cell(:,2)
      write(main_output,*) cell(:,3)
      do i = 1,natoms
        write(main_output,'(i5,2x,a4,2x,6f10.5,4i4)') i,atom_name(i),  &
                  xf(i),yf(i),zf(i),xo(i),yo(i),zo(i), &
                  reference_number(i),reference_cell(i,:)
      end do
    end if

    if (lsort) call sort_structure

    return

  end subroutine rectangulate


!===============================================================================

    subroutine reduce
!   =================

!--------------------------------------------------------------------------
!
!     This subroutine moves all atoms within a supercell configuration back
!     into the original unit cell
!
!--------------------------------------------------------------------------

    use structure_data
    use channel_numbers

!   Coordinates come in as xf(:), with natoms as the number of atoms
!   Cartesian coordinates come in as xo(:)
!   Lattice parameters come in as a,b,c,alpha,beta,gamma
!   Cell vectors come in as cell(3,3); we need to be careful as to which index is
!   which.
!   Supercell integers come in as ncell(3)
!   Trivially, new unit cell parameters are a/ncell(1) etc; I think it is fine
!   to define these as new because we will need a new structure defined in order
!   to write out a final CIF file.
!
!   If you scale down the unit cell, you then need to scale up the fractional
!   coordinates. Thus you should have statements such as xf(i) = xf(i) * ncell(1)
!   etc. But now your fraction coordinate is too large, so you need to subtract
!   the part greater than 1, for example xf(i) = xf(i) - int(xf(i)). Does this
!   work if xf(i) is negative?
!
!   For cartesian coordinates, you need to subtract the components of cell(3,3)
!   from the cartesian coordinates xo(i). The temptation on day 1 is not to bother
!   because in reality we are not likely to need these. But we should set this up
!   so that we can do this one day. There are several options, one of which is to
!   use the fractional coordinates subtraction to work out how many unit cell
!   vectors to subtract.

    use structure_data

    implicit none

    double precision :: al,be,ga,s,volume

    if ((ncell(1)==1).and.(ncell(2)==1).and.(ncell(3)==1)) then
      write(6,'(a)') 'No supercell has been given for structure reduce and therefore no reduction will be made'
      return
    end if
    
    xf = xf*ncell(1) + 1.5d0  ;  xf = xf - aint(xf) - 0.5d0
    yf = yf*ncell(2) + 1.5d0  ;  yf = yf - aint(yf) - 0.5d0
    zf = zf*ncell(3) + 1.5d0  ;  zf = zf - aint(zf) - 0.5d0

    a = a/ncell(1)  ;  b = b/ncell(2)  ;  c = c/ncell(3)
    al = alpha/57.295779512d0
    be = beta/57.295779512d0
    ga = gamma/57.295779512d0
    s = 0.5*(al+be+ga)
!    volume = 2.0*a*b*c*sqrt(sin(s)*sin(s-al)*sin(s-be)*sin(s-ga))
    volume = a*b*c*sqrt(1.0d0-cos(2.0d0*al)-cos(2.0d0*be)-cos(2.0d0*ga)) + 2.0d0*(cos(al)*cos(be)*cos(ga))
    cell(1,1) = volume/(sin(al)*b*c)
    cell(2,1) = a*cos(ga)*sin(al)
    cell(3,1) = a*cos(be)
    cell(1,2) = 0.0d0
    cell(2,2) = b*sin(al)
    cell(3,2) = b*cos(al)
    cell(1,3) = 0.0d0
    cell(2,3) = 0.0d0
    cell(3,3) = c

    xo = cell(1,1)*xf + cell(1,2)*yf + cell(1,3)*zf
    yo = cell(2,1)*xf + cell(2,2)*yf + cell(2,3)*zf
    zo = cell(3,1)*xf + cell(3,2)*yf + cell(3,3)*zf
    
    return

  end subroutine reduce


!===============================================================================

    subroutine nanoparticle
!   =======================

!--------------------------------------------------------------------------
!
!     This subroutine creates a nanoparticle of a chosen shape from the expanded
!     configuration
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use channel_numbers

    implicit none

    integer :: i,j,newatoms,iel
    integer, allocatable :: atom_typet(:),atom_labelt(:),rnumt(:),rcellt(:,:)
    double precision,allocatable :: xft(:),yft(:),zft(:),xot(:),yot(:),zot(:)
    double precision :: rcentre(3),rsq,r,rx,ry,rz
    character(len=10), allocatable :: catom_labelt(:)
    character(len=4), allocatable :: atom_namet(:)
    logical, allocatable :: lwithin(:)

    allocate(xft(natoms),yft(natoms),zft(natoms))
    allocate(xot(natoms),yot(natoms),zot(natoms))
    allocate(atom_namet(natoms))
    allocate(atom_typet(natoms))
!    if (luselabels) then
!      allocate(catom_labelt(natoms))
!    else
!      allocate(atom_labelt(natoms))
!    end if

 
    allocate(atom_labelt(natoms))
    allocate(rnumt(natoms),rcellt(natoms,3))
    
    if(lrotate) then
      xft = xf(1:natoms)
      yft = yf(1:natoms)
      zft = zf(1:natoms)
      xot = xrot(1:natoms)
      yot = yrot(1:natoms)
      zot = zrot(1:natoms)
    else
      xft = xf(1:natoms)
      yft = yf(1:natoms)
      zft = zf(1:natoms)
      xot = xo(1:natoms)
      yot = yo(1:natoms)
      zot = zo(1:natoms)
    end if
    
    atom_namet = atom_name(1:natoms)
    atom_typet = atom_type(1:natoms)
    
!    if (luselabels) then
!      catom_labelt = catom_label(1:natoms)
!    else
!      atom_labelt = atom_label(1:natoms)
!    end if
    atom_labelt = atom_label(1:natoms)
    rnumt = reference_number(1:natoms)
    rcellt = reference_cell(1:natoms,:)
    if (allocated(xf)) deallocate(xf,yf,zf)
    if (allocated(xo)) deallocate(xo,yo,zo)
    if (allocated(atom_name)) deallocate(atom_name)
    if (allocated(atom_type)) deallocate(atom_type)
    if (allocated(atom_label)) deallocate(atom_label)
    if (allocated(catom_label)) deallocate(catom_label)
    if (allocated(reference_number)) deallocate(reference_number)
    if (allocated(reference_cell)) deallocate(reference_cell)

!    rcentre(1) = (cell(1,1) + cell(1,2) + cell(1,3))/2.0d0
!    rcentre(2) = (cell(2,1) + cell(2,2) + cell(2,3))/2.0d0
!    rcentre(3) = (cell(3,1) + cell(3,2) + cell(3,3))/2.0d0
    rcentre(1) = (cell(1,1) + cell(2,1) + cell(3,1))/2.0d0
    rcentre(2) = (cell(1,2) + cell(2,2) + cell(3,2))/2.0d0
    rcentre(3) = (cell(1,3) + cell(2,3) + cell(3,3))/2.0d0
    allocate(lwithin(natoms))
    lwithin = .false.
    newatoms = 0

    do i = 1,natoms
      select case(trim(cshape))
        case('sphere')
         if(lhollow) then
          rx = xot(i) - rcentre(1)
          ry = yot(i) - rcentre(2)
          rz = zot(i) - rcentre(3)
          rsq = rx**2 + ry**2 + rz**2
          r = sqrt(rsq)
          if ((rhollow(1)/2.0d0<=r).and.(r<=rnano(1)/2.0d0)) then
          lwithin(i) = .true.
          newatoms = newatoms + 1
          end if
         else
      rx = xot(i) - rcentre(1)
      ry = yot(i) - rcentre(2)
      rz = zot(i) - rcentre(3)
      rsq = rx**2 + ry**2 + rz**2
      r = sqrt(rsq)
          if (r<=rnano(1)/2.0d0) then
          lwithin(i) = .true.
          newatoms = newatoms + 1
          end if
         end if
        case('cube')
         if(lhollow) then
          rx = xot(i) - rcentre(1)
          ry = yot(i) - rcentre(2)
          rz = zot(i) - rcentre(3)
          if ((abs(rx)<=rnano(1)/2.0d0).and.(abs(ry)<=rnano(2)/2.0d0).and.(abs(rz)<=rnano(3)/2.0d0).and.&
          ((rhollow(1)/2.0d0<=abs(rx)).or.(rhollow(2)/2.0d0<=abs(ry)).or.(rhollow(3)/2.0d0<=abs(rz))))then
          lwithin(i) = .true.
          newatoms = newatoms + 1
          end if 
         else
          rx = xot(i) - rcentre(1)
          ry = yot(i) - rcentre(2)
          rz = zot(i) - rcentre(3)
          if ((abs(rx)<=rnano(1)/2.0d0).and.(abs(ry)<=rnano(2)/2.0d0).and.(abs(rz)<=rnano(3)/2.0d0)) then
          lwithin(i) = .true.
          newatoms = newatoms + 1
          end if
         end if
        case('cylinder')
         if(lhollow) then
          rx = xot(i) - rcentre(1)
          ry = yot(i) - rcentre(2)
          rz = zot(i) - rcentre(3)
          rsq = rx**2 + ry**2
          r = sqrt(rsq)
          if ((abs(rz)<=rnano(3)/2.0d0).and.(r<=rnano(1)/2.0d0).and.&
          (abs(rz)>=rhollow(3)/2.0d0).and.(r>=rhollow(1)/2.0d0)) then
          lwithin(i) = .true.
          newatoms = newatoms + 1
          end if
         else
          rx = xot(i) - rcentre(1)
          ry = yot(i) - rcentre(2)
          rz = zot(i) - rcentre(3)
          rsq = rx**2 + ry**2
          r = sqrt(rsq)
          if ((abs(rz)<=rnano(3)/2.0d0).and.(r<=rnano(1)/2.0d0)) then
          lwithin(i) = .true.
          newatoms = newatoms + 1
          end if
         end if
        case('rectangle')
         if(lhollow) then
          rx = xot(i) - rcentre(1)
          ry = yot(i) - rcentre(2)
          rz = zot(i) - rcentre(3)
          if ((abs(rx)<=rnano(1)/2.0d0).and.(abs(ry)<=rnano(2)/2.0d0).and.(abs(rz)<=rnano(3)/2.0d0).and.&
          ((rhollow(1)/2.0d0<=abs(rx)).or.(rhollow(2)/2.0d0<=abs(ry)).or.(rhollow(3)/2.0d0<=abs(rz)))) then
          lwithin(i) = .true.
          newatoms = newatoms + 1
          end if
         else
          rx = xot(i) - rcentre(1)
          ry = yot(i) - rcentre(2)
          rz = zot(i) - rcentre(3)
          if ((abs(rx)<=rnano(1)/2.0d0).and.(abs(ry)<=rnano(2)/2.0d0).and.(abs(rz)<=rnano(3)/2.0d0)) then
          lwithin(i) = .true.
          newatoms = newatoms + 1
          end if
         end if
        case('ellipsoid')
         if(lhollow) then
          rx = xot(i) - rcentre(1)
          ry = yot(i) - rcentre(2)
          rz = zot(i) - rcentre(3)
          if (((rx**2/(rnano(1)/2.0d0)**2)+(ry**2/(rnano(2)/2.0d0)**2)+(rz**2/(rnano(3)/2.0d0)**2)<=1).and.&
          ((rx**2/(rhollow(1)/2.0d0)**2)+(ry**2/(rhollow(2)/2.0d0)**2)+(rz**2/(rhollow(3)/2.0d0)**2)>=1)) then
        lwithin(i) = .true.
        newatoms = newatoms + 1
      end if
         else
          rx = xot(i) - rcentre(1)
          ry = yot(i) - rcentre(2)
          rz = zot(i) - rcentre(3)
          if ((rx**2/(rnano(1)/2.0d0)**2)+(ry**2/(rnano(2)/2.0d0)**2)+(rz**2/(rnano(3)/2.0d0)**2)<=1) then
          lwithin(i) = .true.
          newatoms = newatoms + 1
          end if
         end if
      end select
    end do
    write(6,'(a,i0)') 'Number of atoms within nanoparticle is ',newatoms

    allocate(xo(newatoms),yo(newatoms),zo(newatoms))
    allocate(xf(newatoms),yf(newatoms),zf(newatoms))
    allocate(atom_name(newatoms))
    allocate(atom_type(newatoms))
!    if (luselabels) then
!      allocate(catom_label(newatoms))
!    else
!      allocate(atom_label(newatoms))
!    end if
    allocate(atom_label(newatoms))
    allocate(reference_number(newatoms))
    allocate(reference_cell(newatoms,3))
    cell = 0.0d0
    
    if (lnanobox_list) then
       ! >>>>>>>>>>>>> Yuanpeng -> Add in the functionality for user to input their own >>>>>>
       ! box size for the nanoparticle. >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       write(*,*) "Warning> If not sure what the box size should be, just leave '-nanobox' out."
       write(*,*) "Warning> Then the default box size multipliers [4 4 4] will be used."
       cell(1,1) = nanobox(1)*rnano(1) ; cell(2,2) = nanobox(2)*rnano(2) ; cell(3,3) = nanobox(3)*rnano(3)
       a = nanobox(1)*rnano(1) ; b = nanobox(2)*rnano(2) ; c = nanobox(3)*rnano(3)
    else
       ! >>>>>>>>>>>>> Yuanpeng -> Here the cell size is set to be four times the dimension of the >>>>
       ! generated nanoparticle, and the reason is to leave enough void space to prevent atoms from >>>
       ! seeing any mirrors of any atoms in the system. In such a way, we are creating a really >>>>>>>
       ! isolated nanoparticle. >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       cell(1,1) = 4.0d0*rnano(1) ; cell(2,2) = 4.0d0*rnano(2) ; cell(3,3) = 4.0d0*rnano(3) !Ask Martin
       a = 4.0d0*rnano(1) ; b = 4.0d0*rnano(2) ; c = 4.0d0*rnano(3)
    end if
    alpha = 90.0d0 ; beta = 90.0d0 ; gamma = 90.0d0
    ! >>>>>>>>>>>> Yuanpeng -> Correct the number density for the big box. >>>
    density = newatoms/(cell(1,1)*cell(2,2)*cell(3,3))
    j = 0
    do i = 1,natoms
      if (lwithin(i)) then
        j = j + 1
        xo(j) = xot(i) ; yo(j) = yot(i) ; zo(j) = zot(i)
        atom_name(j) = atom_namet(i)
        atom_type(j) = atom_typet(i)
!        if (luselabels) then
!          catom_label(j) = catom_labelt(i)
!        else
!          atom_label(j) = atom_labelt(i)
!        end if
        atom_label(j) = atom_labelt(i)
        reference_number(j) = rnumt(i)
        reference_cell(j,:) = rcellt(i,:)
      end if
    end do
    xo = xo - sum(xo)/dble(newatoms) + cell(1,1)/2
    yo = yo - sum(yo)/dble(newatoms) + cell(2,2)/2
    zo = zo - sum(zo)/dble(newatoms) + cell(3,3)/2
    xf = xo/cell(1,1) ; yf = yo/cell(2,2) ; zf = zo/cell(3,3)


    if (allocated(xft)) deallocate(xft,yft,zft)
    if (allocated(xot)) deallocate(xot,yot,zot)
    if (allocated(atom_namet)) deallocate(atom_namet)
    if (allocated(atom_typet)) deallocate(atom_typet)
    if (allocated(atom_labelt)) deallocate(atom_labelt)
    if (allocated(catom_labelt)) deallocate(catom_labelt)
    if (allocated(rnumt)) deallocate(rnumt)
    if (allocated(rcellt)) deallocate(rcellt)

    if (ldiag) then
      write(main_output,*)
      write(main_output,'(a)') 'Output from nanoparticle'
      write(main_output,'(a)') ' === New structure ==='
      write(main_output,'(a)') ' New lattice vectors'
      write(main_output,*) cell(:,1)
      write(main_output,*) cell(:,2)
      write(main_output,*) cell(:,3)
      do i = 1,newatoms
         write(main_output,'(i5,2x,a4,2x,6f10.5,4i4)') i,atom_name(i),  &
                  xf(i),yf(i),zf(i),xo(i),yo(i),zo(i), &
                  reference_number(i),reference_cell(i,:)
      end do
    end if

    natoms = newatoms

    call assign_elements
    ntypes = 0
    if (allocated(numoftype)) then
      do iel = -1,nelements
        if (n_elements(iel)>0) then
          ntypes = ntypes + 1
          numoftype(ntypes) = n_elements(iel)
        end if
      end do
    end if

    if (lsort) call sort_structure

    return

  end subroutine nanoparticle

!===============================================================================

  subroutine scale_nanoparticle
! =============================

!--------------------------------------------------------------------------
!
!     This subroutine scales the size of the inorganic code of a nanoparticle 
!     and shifts the molecular ligands radially outwards
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use channel_numbers

    implicit none

    integer :: i,j,newatoms,iel
    integer, allocatable :: atom_typet(:),atom_labelt(:),rnumt(:),rcellt(:,:)
    double precision,allocatable :: xft(:),yft(:),zft(:),xot(:),yot(:),zot(:)
    double precision :: rcentre(3),rsq,r,rx,ry,rz
    character(len=10), allocatable :: catom_labelt(:)
    character(len=4), allocatable :: atom_namet(:)
    logical, allocatable :: lwithin(:)

    allocate(xft(natoms),yft(natoms),zft(natoms))
    allocate(xot(natoms),yot(natoms),zot(natoms))
    allocate(atom_namet(natoms))
    allocate(atom_typet(natoms))
    allocate(atom_labelt(natoms))
    allocate(rnumt(natoms),rcellt(natoms,3))

! First task is to identify the centre of the nanoparticle core

! Second task is to scale the fractional and orthogonal coordinates of the core atoms

! Third task is to identify the centre of each ligand, then shift this centre outwards,
! and then recalculate the positions of the ligands

    return

  end subroutine scale_nanoparticle


!===============================================================================

  subroutine create_vacancies
! ===========================

!--------------------------------------------------------------------------
!
!     This subroutine generates vacancies from list supplied by the user.
!     Here we are looking for strings of the form "Mg2 0.5 : O3 0.4", or
!     "Na 0.5 : O 0.3", where the atom label has an optional label number,
!     and the second number is the fraction of vacant sites. Vacant atoms
!     will be replaced by the Va vacancy symbol.
!
!--------------------------------------------------------------------------

    use arguments
    use annotations
    use structure_data
    use channel_numbers
    
    implicit none
    
    double precision :: fraction,xrand
    integer :: i,n,isite,nsites,ilabel,ie,na,nselected,iatom,nvacancies, &
               mvacancies
    integer, allocatable :: final_list(:)
    character(len=120) :: cvl
    character(len=2) :: ca
    character(len=10) :: clabel
    logical, allocatable :: lselect(:),lvacant(:)
    
    allocate(lselect(natoms),lvacant(natoms))
    
    call remove_spaces(cvacancy_list)
    cvacancy_list = adjustl(cvacancy_list)
    
    nsites = 0
    cvl = cvacancy_list
    n = 1
    do while (n>0)
      n = index(cvl,':')
      cvl = cvl(n+1:)
      nsites = nsites + 1
    end do

    cvl = cvacancy_list
    site_loop: do isite = 1,nsites
      lselect = .false. ; lvacant = .false.
      n = index(cvl,':')
      if (n>0) then
        if (cvl(n-1:n-1)/=' ') then
          cvl(n+1:) = cvl(n:)
          cvl(n:n) = ' '
        end if
      end if
      ca = '  '
      if ((ichar(cvl(2:2))>=ichar('a')).and.(ichar(cvl(2:2))<=ichar('z'))) then
        ca = cvl(1:2)
        cvl = cvl(3:)
      else
        ca(1:1) = cvl(1:1)
        cvl = cvl(2:)
      end if 
      ie = nelements + 1
      element_search: do i = -1,nelements
        if (ca==element(i)) then
          ie = i
          exit element_search
        end if
      end do element_search
      if (ie==(nelements+1)) cycle site_loop
      if (cvl(1:1)==' ') then
        ilabel = 0
        read(cvl,*) fraction
      else
!        if (luselabels) then
!          read(cvl,*) ilabel,fraction
!        else
!          read(cvl,*) clabel,fraction
!        end if
        read(cvl,*) ilabel,fraction
      end if
      n = index(cvl,':')
      if (n>0) cvl = adjustl(cvl(n+1:))
      na = 0
      do iatom = 1,natoms
        if ((ilabel==0).and.(ie==atom_type(iatom))) lselect(iatom) = .true.
        if ((ilabel/=0).and.(ilabel==atom_label(iatom)).and.(ie==atom_type(iatom))) lselect(iatom) = .true.
!        if ((clabel==catom_label(iatom)).and.(ie==atom_type(iatom))) lselect(iatom) = .true.
      end do
      nselected = count(lselect)
      nvacancies = nint(dble(nselected)*fraction)
      mvacancies = 0
      vacancy_loop: do iatom = 1,natoms
        if (lselect(iatom)) then
          call random_number(xrand)
          if (xrand<=fraction) then
            mvacancies = mvacancies + 1
            lvacant(iatom) = .true.
            lselect(iatom) = .false.
            if (mvacancies==nvacancies) exit vacancy_loop
          end if
        end if
      end do vacancy_loop      
      nselected = count(lselect)
      allocate (final_list(nselected))
      i = 0
      do iatom = 1,natoms
        if (lselect(iatom)) then
          i = i + 1
          final_list(i) = iatom
        end if
      end do
      last_few_vacancies_loop: do while (mvacancies<nvacancies)
        call random_number(xrand) ; xrand = xrand*dble(nselected) ; i = int(xrand)+1
        if (i>nselected) i = nselected
        if (lselect(final_list(i))) then
          mvacancies = mvacancies + 1
          lselect(final_list(i)) = .false.
          lvacant(final_list(i)) = .true.
        end if
      end do last_few_vacancies_loop
      deallocate(final_list)
      assign_vacancies: do iatom = 1,natoms
        if (.not.lvacant(iatom)) cycle assign_vacancies
        atom_type(iatom) = 0
        atom_name(iatom) = 'Va'
! This bit needs fixing for luselabels
        if (nsites>1) then
          atom_label(iatom) = isite
        else
          atom_label(iatom) = 0
        end if
      end do assign_vacancies
      n_elements(0) = n_elements(0) + nvacancies
      n_elements(ie) = n_elements(ie) - nvacancies
    end do site_loop
    
    if (nsites>1) n_labels(0) = nsites

    if (ldiag) then
      write(main_output,'(/a)') 'Output from create_vacancies'
      write(main_output,'(a/)') '============================'
      do i = 1,natoms
        write(main_output,'(i4,2x,a4,2x,i3,3f10.4,4i4)') i,atom_name(i),atom_type(i), &
                xf(i),yf(i),zf(i),reference_number(i),reference_cell(i,:)
      end do
      write(main_output,*)
      write(main_output,'(a/)') 'New numbers of atom types'
      do i = -1,nelements
        if (n_elements(i)>0) write(main_output,'(i0,2x,a2,2(2x,i0))') i,element(i),n_elements(i),n_labels(i)
      end do
      write(main_output,*)
    end if

    return
    
  end  subroutine create_vacancies

!===============================================================================

     subroutine split_occupancy
!  ==============================

!--------------------------------------------------------------------------
!
!    This subrutine detects split site occupancy read from cif and tbl file.
!    As a result it will create an argumet for the 'replace_atoms' subroutine
!    and rise flag 'lreplace' to let the program do the replacement.
!    It also checks if the flag 'lreplace' is already used and ask the user
!    for the decision.
!
!--------------------------------------------------------------------------
    
    use structure_data
    use arguments
    
    implicit none
    
    integer :: i,j,k,nreplace,input,ireplace
    integer, allocatable :: replacement_list(:,:)       ! stores atom pairs which has split occupancy
    real, allocatable    :: replacement_list_occ(:)     ! stores occupancy factors for each pair
    real, allocatable    :: replacement_table(:,:)        ! stores occupancy factors for all pair
    logical, allocatable :: replacement_use_rows(:)         ! stores duplicated rows for replacement table
    logical              :: split,doubled,row_equal
    character(len=2)     :: atom_name_1, atom_name_2
    character(len=80)    :: creplacement_list_new
    character(len=80),allocatable    :: replacement_strings(:)
    character(len=80)    :: replacement_string_temp
    real                 :: occ_sum

!    integer :: i1,i2   ! first and second site of an occupancy issue detected

!   Detect it there is a split occupancy of vacancy present
    nreplace = 0
    split = .false.
    creplacement_list_new = ''
    allocate(replacement_table(natoms,natoms))
    allocate(replacement_use_rows(natoms))
    replacement_use_rows = .true.
    replacement_table = 0
    
!   Rewriting generation of the replace command    
    do i = 1, natoms
        do j = 1, natoms
            if ((xf(i)==xf(j)).and.(yf(i)==yf(j)).and.(zf(i)==zf(j))) replacement_table(i,j) = occ(j)
        enddo
    enddo
    
!   Checking duplicates    
     do i = 1, natoms
        do j = i+1, natoms
            row_equal = .false.
            do k = 1, natoms
                row_equal = (row_equal.or.(replacement_table(i,k).ne.replacement_table(j,k)))
            enddo
            replacement_use_rows(j) = replacement_use_rows(j).and.row_equal
        enddo
     enddo

!    Count the numbers of replacemnets     
     nreplace = 0
    do i = 1, natoms
      if (replacement_use_rows(i)) then
          occ_sum = 0
          do j = 1, natoms
              if (replacement_table(i,j) > 0.0) then
                  occ_sum = occ_sum + replacement_table(i,j)
                  nreplace = nreplace + 1
              endif
          enddo
          if (occ_sum==1.0) then
              nreplace = nreplace - 1
              replacement_table(i,i) = 0.0 
          elseif (occ_sum==0.0) then
              replacement_table(i,i) = 1.0
              nreplace = nreplace + 1
          else
              replacement_table(i,i) = 1.0 - occ_sum
          endif
      endif
    enddo  
    
    allocate(replacement_list(nreplace,2))
    allocate(replacement_list_occ(nreplace))
    allocate(replacement_strings(nreplace))
    replacement_list = 0
    replacement_list_occ = 0.0
    
    ireplace = 1
    do i = 1, natoms
      if (replacement_use_rows(i)) then
          occ_sum = 1.0
          do j = 1, natoms
              if (replacement_table(i,j) > 0.0) then
                  replacement_list(ireplace,1) = i
                  if (i==j) then 
                      replacement_list(ireplace,2) = 0
                  else 
                      replacement_list(ireplace,2) = j
                  endif
                  replacement_list_occ(ireplace) = replacement_table(i,j) / occ_sum
                  occ_sum = occ_sum - replacement_list_occ(ireplace)
                  ireplace = ireplace + 1
              endif
          enddo
      endif
    enddo 
    
!   Detect if there are not same atoms with different occupancies which is not supported now    
    do i = 1, natoms
        do j = 1, i-1
            if (atom_name(i) == atom_name(j).and.occ(i) /= occ(j)) then
                write(6,'(a,a,a)') 'Based on the structure file, the site with split occupancy has been detected' 
                write(6,'(a,i3,a,i3)')  'atoms at positions ', j, ' and ', i
                write(6,'(a,a,a)') 'have different occupancies which is not supported at the moment. We need to STOP'
                stop
            endif
        enddo
    enddo
    
    
!   Create flags and argument for -replace command
     if (nreplace > 0) then
        do i = 1, nreplace
         if (replacement_list(i,1)>0) then 
             atom_name_1 = atom_name(replacement_list(i,1)) 
         else
             atom_name_1 = 'Va' 
         end if
         if (replacement_list(i,2)>0) then 
             atom_name_2 = atom_name(replacement_list(i,2)) 
         else 
             atom_name_2 = 'Va' 
         end if            
         write (replacement_string_temp, "(A2, 1x, F10.7, 1x, A2)") atom_name_1, replacement_list_occ(i), atom_name_2
!WS Checking if the replacement_list_new does not contins duplicated 
         doubled = .false.
         do j = 1, i-1
             !WS Extra check for different sites of the same atoms to be detected -> program should then stop
             if ((replacement_strings(j)(1:2) == replacement_string_temp(1:2)).and.(replacement_strings(j)(15:16) & 
                        == replacement_string_temp(15:16)).and. &
             replacement_list_occ(i) /= replacement_list_occ(j)) then
                 write(6,'(a,a,a)') 'Based on the structure file, the site with split occupancy has been detected' 
                 write(6,'(a,a,a)') trim(replacement_strings(j)), ' and ', trim(replacement_string_temp)
                 write(6,'(a,a,a)') 'This kind of split occupancy is not supported at the moment. We need to STOP'
                 stop
             end if
             doubled = doubled.or.replacement_strings(j) == replacement_string_temp
         enddo
         if (.not.doubled) then
             replacement_strings(i) = replacement_string_temp
         else
             replacement_list_occ(i) = 0.0
         end if
        end do
        
        write (creplacement_list_new, "(A)") replacement_strings(1)
        do j=2, nreplace
            if (replacement_list_occ(j) /= 0.0) then
                write (creplacement_list_new, "(A, 1x, A1, 1x, A)") trim(creplacement_list_new),':', trim(replacement_strings(j))
            end if
        enddo
               
        
        if (lreplace) then
           write(6,'(a)') 'The atomic positions obtained from the input file show split site occupancy &
                                or only partial occupancy. The program will proceed with the partial occupancy.'
           write(6,'(a)') 'On the other hand you intend to use -replace argument. Please decide which shoud be used'
           write(6,'(a,a,a)') '    (from structure file)   (1) -replace[', trim(creplacement_list_new), ']' 
           write(6,'(a,a,a)') '    (from argument list)    (2) -replace[', trim(creplacement_list), ']'
           write(6,'(a,a,a)') 'Please select (1) or (2)    >' 
           read(5,*) input
           if (input==1) then
              creplacement_list = creplacement_list_new
           end if
        else
 !   Write a message for user about the split occupancy recognised
            lreplace = .true.
            creplacement_list = creplacement_list_new
            write(6,'(a)') 'The atomic positions obtained from the input file show split site occupancy or & 
                                only partial occupancy. The program will proceed with the partial occupancy.'
            write(6,'(a)') 'Next time you can also use the following argument:'
            write(6,'(a,a,a)') '-replace[', trim(creplacement_list), ']'                   
        end if
     end if
     
    
    if (allocated(replacement_list)) deallocate(replacement_list)
    if (allocated(replacement_strings)) deallocate(replacement_strings)
    if (allocated(replacement_list_occ)) deallocate(replacement_list_occ)
    if (allocated(replacement_table)) deallocate(replacement_table)
    
    return

    end subroutine split_occupancy
     
!===============================================================================

  subroutine replace_atoms
! ========================

!--------------------------------------------------------------------------
!
!     This subroutine allows you to replace some atoms at random by specified
!     atoms, allowing for site disorder where all sites are fully occupied
!     with fractional occupancy of two different atoms. The argment is given
!     as a string within [...] brackets, where the string has the form
!     "Si 0.5 Al", where we assume that the sum of occupancies is one. It
!     is possible to use labels, eg "Ca1 0.5 Mg". Multiple substitutions are
!     separated by commands, eg "Si 0.125 Al : Mg 0.5 Ca".
!
!--------------------------------------------------------------------------

    use arguments
    use annotations
    use structure_data
    use channel_numbers
    
    implicit none
    
    double precision :: fraction,xrand
    integer :: i,n,isite,nsites,ilabel,ie1,ie2,na,nselected,iatom, &
               nreplacements,mreplacements
    integer, allocatable :: final_list(:)
    character(len=120) :: crl
    character(len=10) :: clabel
    character(len=2) :: ca1,ca2
    logical, allocatable :: lselect(:),lreplacement(:)
    
    allocate(lselect(natoms),lreplacement(natoms))
    
    call remove_spaces(creplacement_list)
    creplacement_list = adjustl(creplacement_list)
    
    nsites = 0
    crl = creplacement_list
    n = 1
    do while (n>0)
      n = index(crl,':')
      crl = crl(n+1:)
      nsites = nsites + 1
    end do

    crl = creplacement_list
    site_loop: do isite = 1,nsites
      lselect = .false. ; lreplacement = .false.
      n = index(crl,':')
      if (n>0) then
        if (crl(n-1:n-1)/=' ') then
          crl(n+1:) = crl(n:)
          crl(n:n) = ' '
        end if
      end if
      ca1 = '  ' ; ca2 = '  '
      if ((ichar(crl(2:2))>=ichar('a')).and.(ichar(crl(2:2))<=ichar('z'))) then
        ca1 = crl(1:2)
        crl = crl(3:)
      else
        ca1(1:1) = crl(1:1)
        crl = crl(2:)
      end if 
      ie1 = nelements + 1
      element_search1: do i = -1,nelements
        if (ca1==element(i)) then
          ie1 = i
          exit element_search1
        end if
      end do element_search1
      if (ie1==(nelements+1)) cycle site_loop
      if (crl(1:1)==' ') then
        ilabel = 0
        crl = adjustl(crl)
        read(crl,*) fraction
        n = index(crl,' ') ; crl = adjustl(crl(n:))
      else
        read(crl,*) ilabel,fraction
        n = index(crl,' ') ; crl = adjustl(crl(n:))
        n = index(crl,' ') ; crl = adjustl(crl(n:))
      end if
      if ((ichar(crl(2:2))>=ichar('a')).and.(ichar(crl(2:2))<=ichar('z'))) then
        ca2 = crl(1:2)
        crl = crl(3:)
      else
        ca2(1:1) = crl(1:1)
        crl = crl(2:)
      end if 
      ie2 = nelements + 1
      element_search2: do i = -1,nelements
        if (ca2==element(i)) then
          ie2 = i
          exit element_search2
        end if
      end do element_search2
      if (ie2==(nelements+1)) cycle site_loop
      n = index(crl,':')
      if (n>0) crl = adjustl(crl(n+1:))
      na = 0
      do iatom = 1,natoms
        if ((ilabel==0).and.(ie1==atom_type(iatom))) lselect(iatom) = .true.
        if ((ilabel/=0).and.(ilabel==atom_label(iatom)).and.(ie1==atom_type(iatom))) lselect(iatom) = .true.
!        if ((clabel==catom_label(iatom)).and.(ie1==atom_type(iatom))) lselect(iatom) = .true.
      end do
      nselected = count(lselect)
      nreplacements = nint(dble(nselected)*fraction)
      mreplacements = 0
      replacement_loop: do iatom = 1,natoms
        if (lselect(iatom)) then
          call random_number(xrand)
          if (xrand<=fraction) then
            mreplacements = mreplacements + 1
            lreplacement(iatom) = .true.
            lselect(iatom) = .false.
            if (mreplacements==nreplacements) exit replacement_loop
          end if
        end if
      end do replacement_loop      
      nselected = count(lselect)
!      allocate (final_list(nselected))
      i = 0
      do iatom = 1,natoms
        if (lreplacement(iatom)) then
          i = i + 1
!          final_list(i) = iatom
        end if
      end do
      last_few_replacements_loop: do while (mreplacements<nreplacements)
        call random_number(xrand) ; xrand = xrand*dble(natoms) ; iatom = int(xrand)+1
        if (iatom>natoms) cycle
        if (lreplacement(iatom)) cycle
        if (lselect(iatom)) then
!        if (lselect(final_list(i))) then
          i = i + 1
!          final_list(i) = iatom 
          mreplacements = mreplacements + 1
          lselect(iatom) = .false.
          lreplacement(iatom) = .true.
        end if
      end do last_few_replacements_loop
!      deallocate(final_list)

      assign_replacements: do iatom = 1,natoms
        if (.not.lreplacement(iatom)) cycle assign_replacements
        atom_type(iatom) = ie2
        atom_name(iatom) = ca2
!        if (nsites>1) then
!          atom_label(iatom) = isite
!        else
!          atom_label(iatom) = 0
!        end if
      end do assign_replacements
      n_elements(ie1) = n_elements(ie1) - nreplacements
      n_elements(ie2) = n_elements(ie2) + nreplacements
    end do site_loop
    
!    if (nsites>1) n_labels(0) = nsites

    if (ldiag) then
      write(main_output,'(/a)') 'Output from replace_atoms'
      write(main_output,'(a/)') '========================='
      do i = 1,natoms
        write(main_output,'(i4,2x,a4,2x,i3,3f10.4,4i4)') i,atom_name(i),atom_type(i), &
                xf(i),yf(i),zf(i),reference_number(i),reference_cell(i,:)
      end do
      write(main_output,*)
      write(main_output,'(a/)') 'New numbers of atom types'
      do i = -1,nelements
        if (n_elements(i)>0) write(main_output,'(i0,2x,a2,2(2x,i0))') i,element(i),n_elements(i),n_labels(i)
      end do
      write(main_output,*)
    end if
    
    ntypes = 0

    return
    
    end  subroutine replace_atoms

!===============================================================================

  subroutine shift_origin
! =======================

!--------------------------------------------------------------------------
!
!     This subroutine allows you to shift the origin in the configuration.
!     It is possible to shift along x, y and z axes.
!
!--------------------------------------------------------------------------

    use arguments
    use annotations
    use structure_data
    use channel_numbers
    
    implicit none
    
    integer :: i
    double precision :: xshift,yshift,zshift
    
    xshift = shiftx/a
    yshift = shifty/b
    zshift = shiftz/c
    
    write(6,'(a)') 'Here I am'
    
    do i = 1,natoms
      if (lshiftx) xf(i) = xf(i) + xshift
      if (lshifty) yf(i) = yf(i) + yshift
      if (lshiftz) zf(i) = zf(i) + zshift
    end do

    do i = 1,natoms
      xo(i) = cell(1,1)*xf(i) + cell(1,2)*yf(i) + cell(1,3)*zf(i)
      yo(i) = cell(2,1)*xf(i) + cell(2,2)*yf(i) + cell(2,3)*zf(i)
      zo(i) = cell(3,1)*xf(i) + cell(3,2)*yf(i) + cell(3,3)*zf(i)
    end do

!    if (nsites>1) n_labels(0) = nsites

    if (ldiag) then
      write(main_output,'(/a)') 'Output from shift_origin'
      write(main_output,'(a/)') '========================'
      do i = 1,natoms
        write(main_output,'(i0,2x,a4,2x,i3,3f10.4,4i4)') i,atom_name(i),atom_type(i), &
                xf(i),yf(i),zf(i),reference_number(i),reference_cell(i,:)
      end do
      write(main_output,*)
    end if

    return
    
    end  subroutine shift_origin

!===============================================================================


  subroutine reset_origin
! =======================

!--------------------------------------------------------------------------
!
!     This subroutine resets the origin placing the centre of mass in the
!     centre of the configuration
!
!--------------------------------------------------------------------------

    use arguments
    use annotations
    use structure_data
    use channel_numbers
    
    implicit none
    
    integer :: i
    double precision :: xshift,yshift,zshift
    
    xshift = sum(xf)/natoms - 0.5d0
    yshift = sum(yf)/natoms - 0.5d0
    zshift = sum(zf)/natoms - 0.5d0
    
    xf = xf - xshift
    yf = yf - yshift
    zf = zf - zshift

    do i = 1,natoms
      xo(i) = cell(1,1)*xf(i) + cell(1,2)*yf(i) + cell(1,3)*zf(i)
      yo(i) = cell(2,1)*xf(i) + cell(2,2)*yf(i) + cell(2,3)*zf(i)
      zo(i) = cell(3,1)*xf(i) + cell(3,2)*yf(i) + cell(3,3)*zf(i)
    end do

!    if (nsites>1) n_labels(0) = nsites

    if (ldiag) then
      write(main_output,'(/a)') 'Output from reset_origin'
      write(main_output,'(a/)') '========================'
      do i = 1,natoms
        write(main_output,'(i0,2x,a4,2x,i3,3f10.4,4i4)') i,atom_name(i),atom_type(i), &
                xf(i),yf(i),zf(i),reference_number(i),reference_cell(i,:)
      end do
      write(main_output,*)
    end if

    return
    
    end  subroutine reset_origin

!===============================================================================





  subroutine assign_atoms_to_use
! ==============================

! This assigns atoms to use based on the provided limits

    use structure_data

    latomuse = .false.
    do i = 1,natoms
      if (xf(i)<xyzlimits(1)) cycle
      if (xf(i)>xyzlimits(2)) cycle
      if (yf(i)<xyzlimits(3)) cycle
      if (yf(i)>xyzlimits(4)) cycle
      if (zf(i)<xyzlimits(5)) cycle
      if (zf(i)>xyzlimits(6)) cycle
      latomuse(i) = .true.
    end do

   return
   
   end subroutine assign_atoms_to_use


!===============================================================================

  subroutine statistics
! =====================

! This calculates some basic statistics on the configuration

    use structure_data
    use arguments
    use channel_numbers

    implicit none
    
    integer :: i,j

     if (ntypes==0) then
       do i = -1,nelements
         if (n_elements(i)>0) ntypes = ntypes + 1
       end do
     end if
    if (.not.allocated(numoftype)) then
      allocate(numoftype(ntypes))
      numoftype = 0
      i = 0
      do j = -1,nelements
        if (n_elements(j)>0) then
          i = i + 1
          numoftype(i) = n_elements(j)
        end if
      end do
    end if

    allocate(element_of_type(ntypes),type_number(ntypes),concentration(ntypes))
    concentration = dble(numoftype)/dble(sum(numoftype))

    j = 0
    do i = -1,nelements
      if (n_elements(i)>0) then
        j = j + 1
        element_of_type(j) = element(i)
        type_number(j) = i
      end if
    end do
    
    allocate(ordertype(natoms))
    do i = 1,natoms
      do j = 1,ntypes
        if (atom_type(i)==type_number(j)) then
          ordertype(i) = j
          exit
        end if
      end do
    end do


    if (ldiag) then
      write(main_output,'(/a)') 'Output from statistics'
      write(main_output,'(a/)') '======================'
      write(main_output,'(a,i0,a)') 'There are ',natoms,' atoms in the configuration count'
      write(main_output,'(a,i0,a)') 'There are ',ntypes,' different elements in the configuration'
      write(main_output,'(a,i0,a)') 'There are ',sum(numoftype),' atoms in the configuration from type count'
      do i = 1,ntypes
        write(main_output,'(a,i0,a)') 'The configuration contains ',numoftype(i), &
                                ' atoms of element type '//trim(element_of_type(i))
      end do
      do i = 1,ntypes
        write(main_output,'(a,f4.2)') 'The concentration of element type '//trim(element_of_type(i))//' is ', &
                                      concentration(i)
      end do
      
      write(main_output,*)
      do i = 1,natoms
        write(main_output,'(i0,2x,i0,2x,a,2x,i0)') i,atom_type(i),element(atom_type(i)),ordertype(i)
      end do
    end if

    if (sum(numoftype)/=natoms) then
      write(6,'(a)') 'There is an inconsistency in the count of the number of atoms so we have to stop'
      stop
    end if

  end subroutine statistics


!===============================================================================

  subroutine annotate
! ===================

    use arguments
    use annotations

    implicit none

    integer :: date_time(8),i,n,ierror
    integer, parameter :: imdata=61
    character (len=10) :: cdt1,cdt2,cdt3
    character (len=120) :: buffer,ubuffer,cmetadata
    logical :: lset1, lset2, lset3, lset4
    lm_title = .true. ; lm_owner = .true. ; lm_material = .true.
    lm_comment = .true. ; lm_source = .true.

    lset1 = (lcssr.or.ldlpoly.or.lrmc3.or.lrmc6f.or.lrmc6o)
    lset2 = (lcssr.or.lrmc6f.or.lrmc6o)
    lset3 = (lrmc6f.or.lrmc6o)
    lset4 = (lcssr.or.ldlpoly.or.lrmc6f.or.lrmc6o)

    config_title = ''
    metadata_owner = ''
    metadata_affiliation = ''
    metadata_material = ''
    metadata_phase = ''
    metadata_formula = ''
    metadata_title = ''
    metadata_purpose = ''
    metadata_keywords = ''
    metadata_temperature = ''
    metadata_pressure = ''
    metadata_note = ''
    metadata_comment = ''
    metadata_source = ''
    metadata_date = ''

    if (lannotate) then
    write(6,'(a)') 'Request for metadata for configuration'
    write(6,'(a)') 'Blank replies are acceptable (if foolish); they just cancel the request'
    if (lset4) then
      if (.not.lm_title) then
        write(6,'(a)', advance="no") '... Title for configuration: '
        read(5,'(a)') config_title
      end if
    end if
    if (lset1.and.(cinput(1:3)/='cif')) then
      write(6,'(a)', advance="no") '... Name of material: '
      read(5,'(a)') metadata_material
    end if
    if (lset2) then
      write(6,'(a)', advance="no") '... Name of configuration owner: '
      read(5,'(a)') metadata_owner
    end if
    if (lset2.and.(cinput(1:3)/='cif')) then
      write(6,'(a)', advance="no") '... Source of structure: '
      read(5,'(a)') metadata_source
    end if
    if (lset3) then
      write(6,'(a)', advance="no") '... Additional notes: '
      read(5,'(a)') metadata_comment
    end if
    end if

    if (lmetadata) then
      open(imdata,file=trim(fmetadata),status='old',form='formatted')
      ierror = 0
      read(imdata,'(a)',iostat=ierror) buffer
      metadata_loop: do while(ierror==0)
        n = len_trim(buffer)
        do i = 1,n                          ! this bit converts to uppercase
           if ((ichar(buffer(i:i))>=ichar('a')).and.(ichar(buffer(i:i))<=ichar('z'))) &
                buffer(i:i) = char(ichar(buffer(i:i)) - ichar('a') + ichar('A'))
        end do
        if (index(buffer,'[METADATA]')>0) then
          read(imdata,'(a)',iostat=ierror) buffer
          if (ierror/=0) stop
          do while(len_trim(buffer)>0)
            n = index(buffer,'=')
            ubuffer = trim(buffer(1:n-1))
            cmetadata = adjustl(buffer(n+1:))
            do i = 1,n-1                          ! this bit converts to uppercase
              if ((ichar(ubuffer(i:i))>=ichar('a')).and.(ichar(ubuffer(i:i))<=ichar('z'))) &
                ubuffer(i:i) = char(ichar(ubuffer(i:i)) - ichar('a') + ichar('A'))
            end do
            if (trim(ubuffer)=='INVESTIGATORS') metadata_owner = trim(cmetadata)
            if (trim(ubuffer)=='AFFILIATION') metadata_affiliation = trim(cmetadata)
            if (trim(ubuffer)=='MATERIAL') metadata_material = trim(cmetadata)
            if (trim(ubuffer)=='PHASE') metadata_phase = trim(cmetadata)
            if (trim(ubuffer)=='CHEMICAL FORMULA') metadata_formula = trim(cmetadata)
            if (trim(ubuffer)=='TITLE') metadata_title = trim(cmetadata)
            if (trim(ubuffer)=='PURPOSE') metadata_purpose = trim(cmetadata)
            if (trim(ubuffer)=='KEYWORDS') metadata_keywords = trim(cmetadata)
            if (trim(ubuffer)=='TEMPERATURE') metadata_temperature = trim(cmetadata)
            if (trim(ubuffer)=='PRESSURE') metadata_pressure = trim(cmetadata)
            if (trim(ubuffer)=='NOTE') metadata_note = trim(cmetadata)
            if (trim(ubuffer)=='COMMENT') metadata_comment = trim(cmetadata)
            read(imdata,'(a)',iostat=ierror) buffer
            if (ierror/=0) exit metadata_loop
          end do
          exit metadata_loop
        end if
        read(imdata,'(a)',iostat=ierror) buffer
      end do metadata_loop
      close(imdata)
    end if

    if (len(trim(config_title))==0) then
    config_title = 'No title given' ; lm_title = .false.
    end if
    if (len(trim(metadata_owner))==0) lm_owner = .false.
    if (len(trim(metadata_material))==0) lm_material = .false.
    if (len(trim(metadata_comment))==0) lm_comment = .false.
    if (len(trim(metadata_source))==0) lm_source = .false.

    call date_and_time (cdt1,cdt2,cdt3,date_time)
    metadata_date(1:2) = cdt1(7:8)
    metadata_date(3:3) = '-'
    metadata_date(4:5) = cdt1(5:6)
    metadata_date(6:6) = '-'
    metadata_date(7:10) = cdt1(1:4)
!    if (len_trim(metadata_owner)>0) write(6,'(a)') trim(metadata_owner)
!    if (len_trim(metadata_affiliation)>0) write(6,'(a)') trim(metadata_affiliation)
!    if (len_trim(metadata_material)>0) write(6,'(a)') trim(metadata_material)
!    if (len_trim(metadata_phase)>0) write(6,'(a)') trim(metadata_phase)
!    if (len_trim(metadata_formula)>0) write(6,'(a)') trim(metadata_formula)
!    if (len_trim(metadata_title)>0) write(6,'(a)') trim(metadata_title)
!    if (len_trim(metadata_purpose)>0) write(6,'(a)') trim(metadata_purpose)
!    if (len_trim(metadata_keywords)>0) write(6,'(a)') trim(metadata_keywords)
!    if (len_trim(metadata_temperature)>0) write(6,'(a)') trim(metadata_temperature)
!    if (len_trim(metadata_pressure)>0) write(6,'(a)') trim(metadata_pressure)
!    if (len_trim(metadata_note)>0) write(6,'(a)') trim(metadata_note)
!    if (len_trim(metadata_comment)>0) write(6,'(a)') trim(metadata_comment)

  end subroutine annotate

!===============================================================================

   subroutine remove_errors(text)
!  ==============================

!--------------------------------------------------------------------------
!
!    Replace all bits of text bounded by round brackets
!
!--------------------------------------------------------------------------

    implicit none

    integer :: n,m
    character*(*) :: text

    do while (index(text,'(')>0)
      n = index(text,'(')
      m = index(text,')')
      text = text(1:n-1)//text(m+1:)
    end do

    return

   end subroutine remove_errors

!===============================================================================

   subroutine remove_space(text)
!  =============================

!--------------------------------------------------------------------------
!
!    Replace all spaces in a string
!
!--------------------------------------------------------------------------

    implicit none

    integer :: n
    character*(*) :: text

    do while (index(trim(text),' ')>0)
      n = index(text,' ')
      text = text(1:n-1)//text(n+1:)
    end do

    return

   end subroutine remove_space

!===============================================================================

   subroutine remove_spaces(text)
!  ==============================

!--------------------------------------------------------------------------
!
!    Replace all groups of spaces in a string by a single space
!
!--------------------------------------------------------------------------

    implicit none

    integer :: n
    character*(*) :: text

    do while (index(trim(text),'  ')>0)
      n = index(text,'  ')
      text = text(1:n)//text(n+2:)
    end do

    return

   end subroutine remove_spaces

!===============================================================================
   
      subroutine remove_spaces_and_other(text)
!  ==============================

!--------------------------------------------------------------------------
!
!    Replace all groups of spaces in a string by a single space and also removes all specified characters
!
!--------------------------------------------------------------------------

    implicit none

    integer :: n
    character*(*) :: text

    do while (index(trim(text),'  ')>0)
      n = index(text,'  ')
      text = text(1:n)//text(n+2:)
    end do
    
     do while (index(trim(text),'/')>0)
      n = index(text,'/')
      text = text(1:n-1)//text(n+1:)
    end do   

    return

   end subroutine remove_spaces_and_other

!===============================================================================

      subroutine remove_stuff(text)
!  ==============================

!--------------------------------------------------------------------------
!
!    Replace all groups of spaces in a string by a single space and also removes all specified characters
!
!--------------------------------------------------------------------------

    implicit none

    integer :: n
    character*(*) :: text

    n = index(text,'  ')
    do while (index(trim(text),'  ')>0)
      n = index(text,'  ')
      text = text(1:n)//text(n+2:)
    end do
    text = text(1:n-1)
    return

   end subroutine remove_stuff


!===============================================================================

  subroutine obtain_structures_history
! ====================================

!--------------------------------------------------------------------------
!
!     This subroutine reads in the configurations from a DL_POLY HISTORY
!     file and generates new files
!
!--------------------------------------------------------------------------

   use structure_data
   use arguments
   use channel_numbers

   implicit none

   integer :: i,ilines,ireadflag,iconfig,nconfigurations,nlines,ierror
   double precision :: volume,x,y,z,rvolume,rcell(3,3)
   double precision :: calpha, cbeta, cgamma
   character(len=80) :: cbuffer,title
   character(len=120) :: filetemp
      
   call obtain_number_of_lines_of_file(nlines,filename)
   
   open(ic,file=trim(filename),form='formatted',status='old')
   read(ic,*)
   read(ic,*) ilines,i,natoms
   
   nconfigurations = nint(dble(nlines - 2)/dble(4+2*natoms*(1+ilines)))
   write(6,'(i0,a)') nconfigurations,' configurations have been identified in the HISTORY file'
   filetemp = fileroot
   
   do iconfig = 1,nconfigurations
	   read(ic,*) 
	   do i=1,3 ; read(ic,*) cell(i,:) ; end do
	   volume = cell(1,1)*cell(2,2)*cell(3,3) + cell(1,2)*cell(2,3)*cell(3,1) + &
				cell(1,3)*cell(2,1)*cell(3,2) - cell(1,1)*cell(2,3)*cell(3,2) - &
				cell(1,3)*cell(2,2)*cell(3,1) - cell(1,2)*cell(2,1)*cell(3,3)
	   a = sqrt(cell(1,1)**2 + cell(1,2)**2 + cell(1,3)**2)
	   b = sqrt(cell(2,1)**2 + cell(2,2)**2 + cell(2,3)**2)
	   c = sqrt(cell(3,1)**2 + cell(3,2)**2 + cell(3,3)**2)
	   calpha = (cell(2,1)*cell(3,1) + cell(2,2)*cell(3,2) + cell(2,3)*cell(3,3))/(b*c)
	   cbeta  = (cell(1,1)*cell(3,1) + cell(1,2)*cell(3,2) + cell(1,3)*cell(3,3))/(a*c)
	   cgamma = (cell(1,1)*cell(2,1) + cell(1,2)*cell(2,2) + cell(1,3)*cell(2,3))/(a*b)
	   alpha = acos(calpha)*90.0d0/asin(1.0d0)
	   beta  = acos(cbeta)*90.0d0/asin(1.0d0)
	   gamma = acos(cgamma)*90.0d0/asin(1.0d0)

	   volume = cell(1,1)*cell(2,2)*cell(3,3) + cell(1,2)*cell(2,3)*cell(3,1) + &
				cell(1,3)*cell(2,1)*cell(3,2) - cell(1,1)*cell(2,3)*cell(3,2) - &
				cell(1,3)*cell(2,2)*cell(3,1) - cell(1,2)*cell(2,1)*cell(3,3)

	   rcell(1,1) = cell(2,2)*cell(3,3) - cell(2,3)*cell(3,2) ! correct
	   rcell(1,2) = cell(2,3)*cell(3,1) - cell(2,1)*cell(3,3) ! correct
	   rcell(1,3) = cell(2,1)*cell(3,2) - cell(2,2)*cell(3,1) ! correct
	   rcell(2,1) = cell(3,2)*cell(1,3) - cell(3,3)*cell(1,2) ! correct
	   rcell(2,2) = cell(3,3)*cell(1,1) - cell(3,1)*cell(1,3) ! correct
	   rcell(2,3) = cell(3,1)*cell(1,2) - cell(3,2)*cell(1,1) ! correct
	   rcell(3,1) = cell(1,2)*cell(2,3) - cell(1,3)*cell(2,2) ! correct
	   rcell(3,2) = cell(1,3)*cell(2,1) - cell(1,1)*cell(2,3) ! correct
	   rcell(3,3) = cell(1,1)*cell(2,2) - cell(1,2)*cell(2,1) ! correct

	   rcell = rcell/volume

	   rvolume = rcell(1,1)*rcell(2,2)*rcell(3,3) + rcell(1,2)*rcell(2,3)*rcell(3,1) + &
				 rcell(1,3)*rcell(2,1)*rcell(3,2) - rcell(1,1)*rcell(2,3)*rcell(3,2) - &
				 rcell(1,3)*rcell(2,2)*rcell(3,1) - rcell(1,2)*rcell(2,1)*rcell(3,3)

	   allocate(xf(natoms),yf(natoms),zf(natoms))
	   allocate(xo(natoms),yo(natoms),zo(natoms))
	   allocate(atom_name(natoms))
	   allocate(atom_type(natoms))
	   allocate(atom_label(natoms))

	   density = dble(natoms)/volume

       ierror = 0
	   do i = 1,natoms
		 read(ic,'(a)') cbuffer
		 cbuffer = adjustl(cbuffer)
		 atom_name(i) = cbuffer(1:2)
		 read(ic,*,iostat=ierror) x,y,z
		 if (ierror/=0) then
		   write(6,*) iconfig,i
		   stop
		 end if
		 xo(i) = x ; yo(i) = y ; zo(i) = z
		 xf(i) = rcell(1,1)*x + rcell(1,2)*y + rcell(1,3)*z
		 yf(i) = rcell(2,1)*x + rcell(2,2)*y + rcell(2,3)*z
		 zf(i) = rcell(3,1)*x + rcell(3,2)*y + rcell(3,3)*z
		 if (xf(i)<0) xf(i) = xf(i) + 1.0d0
		 if (yf(i)<0) yf(i) = yf(i) + 1.0d0
		 if (zf(i)<0) zf(i) = zf(i) + 1.0d0
		 if (xf(i)>1) xf(i) = xf(i) - 1.0d0
		 if (yf(i)>1) yf(i) = yf(i) - 1.0d0
		 if (zf(i)>1) zf(i) = zf(i) - 1.0d0
		 if (ilines>=1) read(ic,*)
		 if (ilines==2) read(ic,*)
	   end do

	   call assign_elements
	   
	  write(fileroot,'(i0)') iconfig
	  fileroot = trim(filename)//'_'//adjustl(fileroot)
	   
      if (lcssr) call write_cssr
      if (lcmtx) call write_cmtx
      if (lcif) call write_cif
      if (latomeye) call write_atomeye
      if (lgasp) call write_gasp
      if (lrmc3) call write_rmc3
      if (lrmc6f) call write_rmc6f
      if (lrmc6o) call write_rmc6o
      if (lcrystal) call write_crystal
	  deallocate(xf,yf,zf)
	  deallocate(xo,yo,zo)
	  deallocate(atom_name)
	  deallocate(atom_type)
	  deallocate(atom_label)

   end do

   return

  end subroutine obtain_structures_history
  
  
!===============================================================================


  subroutine obtain_number_of_lines_of_file(nlines,filename)
! ================================================================

!--------------------------------------------------------------------------
!
!     This subroutine obtains the number of lines in a file
!
!--------------------------------------------------------------------------

  use channel_numbers

  integer :: nlines
  character(len=*) :: filename
  character(len=132) :: buffer
  logical :: lexist

  inquire(file=trim(filename),exist=lexist)
  if (.not.lexist) then
    write(6,'(a)') 'File entitled '//trim(filename)//' not found'
    stop
  end if

  buffer = 'wc -l '//trim(filename)//' > .nlines'
! Needs to check especially on non-Unix based OS.
  call execute_command_line(trim(buffer))
  open(istuff,file='.nlines',form='formatted',status='old')
  read(istuff,*) nlines
  call execute_command_line('rm .nlines')
  
  return
  
  end subroutine obtain_number_of_lines_of_file


!===============================================================================

  subroutine write_cif
! ====================

!--------------------------------------------------------------------------
!
!     This subroutine writes the configuration in cif format, suitable
!     to be read by CrystalMaker or another crystal viewer program
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use annotations
    use version_data
    use channel_numbers

    implicit none

    integer :: i,n
    logical :: lokay
    character(len=132) :: fname

    n = index(fileroot,'.')
    if (n>0) then
      fname = trim(wfolder)//fileroot(1:n)//'cif'
    else
      fname = trim(wfolder)//trim(fileroot)//'.cif'
    end if
    inquire(file=trim(fname),exist=lokay)
    if (lokay) then
      n = index(fname,'.',back=.true.)
      fname = 'new_'//trim(fname)
    end if
    open(icif,file=trim(fname),form='formatted',status='unknown')

    write(icif,'(a)') 'data_'//fileroot(1:n-1)
    write(icif,'(a)') '_audit_creation_method         ''generated by data2config v' &
                      //trim(version_number(nversions))//''''

    write(icif,*)
    write(icif,'(a,16x,f10.5)') '_cell_length_a',a
    write(icif,'(a,16x,f10.5)') '_cell_length_b',b
    write(icif,'(a,16x,f10.5)') '_cell_length_c',c
    write(icif,'(a,13x,f10.5)') '_cell_angle_alpha',alpha
    write(icif,'(a,14x,f10.5)') '_cell_angle_beta',beta
    write(icif,'(a,13x,f10.5)') '_cell_angle_gamma',gamma

    write(icif,*)
    write(icif,'(a)') 'loop_'
    write(icif,'(a)') '_symmetry_equiv_pos_as_xyz'
    write(icif,'(a)') '''+x,+y,+z'''

    write(icif,*)
    write(icif,'(a)') 'loop_'
    write(icif,'(a)') '_atom_site_type_symbol'
    write(icif,'(a)') '_atom_site_label'
    write(icif,'(a)') '_atom_site_fract_x'
    write(icif,'(a)') '_atom_site_fract_y'
    write(icif,'(a)') '_atom_site_fract_z'
    do i = 1,natoms
      write(icif,'(a2,4x,a,i0,4x,3f10.4)') element(atom_type(i)),trim(element(atom_type(i))),i,xf(i),yf(i),zf(i)
    end do
    close(icif)

    return

  end subroutine write_cif

!===============================================================================

  subroutine write_cssr
! =====================

!--------------------------------------------------------------------------
!
!     This subroutine writes the configuration in the cssr format, suitable
!     to be read by CrystalMaker
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use annotations
    use channel_numbers

    implicit none

    integer :: i,n
    character(len=132) :: text

    n = index(fileroot,'.')
    if (n>0) then
      open(icssr,file=trim(wfolder)//fileroot(1:n)//'cssr',form='formatted',status='unknown')
    else
      open(icssr,file=trim(wfolder)//trim(fileroot)//'.cssr',form='formatted',status='unknown')
    endif

    write(icssr,'(38x,3f8.3)') a,b,c
    write(icssr,'(22x,3f8.3,a)') alpha,beta,gamma,'    SPGR =1   P 1'

    text = ''
    if (lm_title) text = adjustl(trim(config_title))
    if (lm_material) text = trim(text)//'; material = '//adjustl(trim(metadata_material))
    write(icssr,'(i6,i2,a)') natoms,0,' '//trim(text)
    text = 'Date = '//adjustl(trim(metadata_date))
    if (lm_owner) text = trim(text)//'; Owner = '//adjustl(trim(metadata_owner))
    if (lm_source) text = trim(text)//'; Source = '//adjustl(trim(metadata_source))
    write(icssr,'(a)') text
    do i = 1,natoms
    write(icssr,'(i6,1x,a2,4x,3f10.4)') i,element(atom_type(i)),xf(i),yf(i),zf(i)
    end do
    close(icssr)

    return

  end subroutine write_cssr


!===============================================================================

  subroutine write_cmtx
! =====================

!--------------------------------------------------------------------------
!
!     This subroutine writes the configuration in the cmtx format, suitable
!     to be read by CrystalMaker
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use annotations
    use channel_numbers

    implicit none

    integer :: i,n,style
    character(len=32) :: atom_label_local
    logical :: lmcfg
    double precision :: mspin(3),rcell(3,3),volume,vlength,vlengthmax,suvwmax,s
    double precision, allocatable :: suvw(:,:)
    
    vlengthmax = 2.0
    style = 2

   volume = cell(1,1)*cell(2,2)*cell(3,3) + cell(1,2)*cell(2,3)*cell(3,1) + &
            cell(1,3)*cell(2,1)*cell(3,2) - cell(1,1)*cell(2,3)*cell(3,2) - &
            cell(1,3)*cell(2,2)*cell(3,1) - cell(1,2)*cell(2,1)*cell(3,3)

   rcell(1,1) = cell(2,2)*cell(3,3) - cell(2,3)*cell(3,2) ! correct
   rcell(1,2) = cell(2,3)*cell(3,1) - cell(2,1)*cell(3,3) ! correct
   rcell(1,3) = cell(2,1)*cell(3,2) - cell(2,2)*cell(3,1) ! correct
   rcell(2,1) = cell(3,2)*cell(1,3) - cell(3,3)*cell(1,2) ! correct
   rcell(2,2) = cell(3,3)*cell(1,1) - cell(3,1)*cell(1,3) ! correct
   rcell(2,3) = cell(3,1)*cell(1,2) - cell(3,2)*cell(1,1) ! correct
   rcell(3,1) = cell(1,2)*cell(2,3) - cell(1,3)*cell(2,2) ! correct
   rcell(3,2) = cell(1,3)*cell(2,1) - cell(1,1)*cell(2,3) ! correct
   rcell(3,3) = cell(1,1)*cell(2,2) - cell(1,2)*cell(2,1) ! correct

   rcell = rcell/volume

    n = index(fileroot,'.')
    if (n>0) then
      open(icmtx,file=trim(wfolder)//fileroot(1:n)//'cmtx',form='formatted',status='unknown')
      inquire(file=trim(wfolder)//fileroot(1:n)//'cfg',exist=lmcfg)
    else
      open(icmtx,file=trim(wfolder)//trim(fileroot)//'.cmtx',form='formatted',status='unknown')
      inquire(file=trim(wfolder)//fileroot(1:n)//'.cfg',exist=lmcfg)
    endif

    write(icmtx,'(a)') 'TITL  New Crystal'
    write(icmtx,'(a,2x,6f12.5/)') 'CELL',a,b,c,alpha,beta,gamma

    write(icmtx,'(a)') '! Lattice type'
    write(icmtx,'(a/)') 'LATC  P'

    write(icmtx,'(a)') '! Spacegroup symbol:  P 1'
    write(icmtx,'(a)') '! Total of 1 general equivalent positions (excl. lattice operators)'
    write(icmtx,'(a/)') 'SYMM     +x       +y       +z  '

    write(icmtx,'(a)') '! Plot range (min/max along x, y and z axes)'
    write(icmtx,'(a/)') 'XYZR  0.000000 1.000000 0.000000 1.000000 0.000000 1.000000'

    write(icmtx,'(a)') '! Model type'
    write(icmtx,'(a/)') 'MODL  1'

    write(icmtx,'(a)') '! Background colour'
    write(icmtx,'(a/)') 'BKCL  1.000 1.000 1.000'

    write(icmtx,'(a)') '! Bond colour'
    write(icmtx,'(a/)') 'BNCL  0.6867 0.6867 0.6867'

    write(icmtx,'(a)') '! Relative bond radii for ball-and-stick and stick plots'
    write(icmtx,'(a/)') 'BRAD    0.2500   0.5000'

    write(icmtx,'(a)') '! Illumination direction'
    write(icmtx,'(a/)') 'LTDN  -0.500000 0.500000 0.707100'

    write(icmtx,'(a)') '! Atom, bond, polyhedral and surface rendering coefficients'
    write(icmtx,'(a)') 'AREN  0.250000 0.750000 1.000000 30'
    write(icmtx,'(a)') 'BREN  0.500000 0.900000 1.000000 10'
    write(icmtx,'(a)') 'PREN  0.500000 0.900000 0.000000 1'
    write(icmtx,'(a/)') 'SREN  0.500000 1.000000 0.000000 1'

    write(icmtx,'(a)') '! Unit cell visibility (1=true; 0 = false)'
    write(icmtx,'(a/)') 'SHCL  1'

!    write(icmtx,'(a)') '! Orientation matrix: 11 12 13; 21 22 23; 31 32 33'
!    write(icmtx,'(a/)') 'OMAT  4.000 0.000 -0.000  -0.000 -5.000 0.000  -0.000 -0.000 -6.000'


    write(icmtx,'(a)') '! Asymmetric unit'
    write(icmtx,'(a)') 'ATOM'
    n = 0
    do i = 1,natoms
      if (.not.latomuse(i)) cycle
      n = n + 1
      write(atom_label_local,'(i0)') n
      atom_label_local = trim(element(atom_type(i)))//adjustl(atom_label_local)
      write(icmtx,'(a2,2x,a10,3f12.6)') element(atom_type(i)),trim(atom_label_local),xf(i),yf(i),zf(i)
    end do
    
    if (lmcfg) then
      write(icmtx,'(/a)') '! Atom Vector: apply vectors on a per-ATOM basis'
      write(icmtx,'(a)') 'AVEC'
      write(icmtx,'(a)') '! <atom label>	  <xFrac>  <yFrac>  <zFrac>	<vecU> <vecV> <vecW> <length>  <red> <green> <blue>    <style>'
      n = index(fileroot,'.')
      if (n>0) then
        open(imcfg,file=trim(wfolder)//fileroot(1:n)//'cfg',form='formatted',status='unknown')
      else
        open(imcfg,file=trim(wfolder)//trim(fileroot)//'.cfg',form='formatted',status='unknown')
      endif
      read(imcfg,*)
      read(imcfg,*) nspins
      allocate(suvw(nspins,3))
      do i = 1,nspins
        read(imcfg,*) mspin
        suvw(i,1) = rcell(1,1)*mspin(1) + rcell(1,2)*mspin(2) + rcell(1,3)*mspin(3)
        suvw(i,2) = rcell(2,1)*mspin(1) + rcell(2,2)*mspin(2) + rcell(2,3)*mspin(3)
        suvw(i,3) = rcell(3,1)*mspin(1) + rcell(3,2)*mspin(2) + rcell(3,3)*mspin(3)
      end do

      suvwmax = 0.0d0
      do i = 1,nspins
        s = sqrt(suvw(i,1)**2 + suvw(i,2)**2 + suvw(i,3)**2)
        if (s > suvwmax) suvwmax = s
      end do

      do i = 1,nspins
        write(atom_label_local,'(i0)') i
        atom_label_local = trim(element(atom_type(i)))//adjustl(atom_label_local)
        s = sqrt(suvw(i,1)**2 + suvw(i,2)**2 + suvw(i,3)**2)
        vlength = vlengthmax*s/suvwmax
        write(icmtx,'(a,2x,6f12.6,f6.2,3f6.3,1x,i0)') trim(atom_label_local), &
                     xf(i),yf(i),zf(i),suvw(i,:),vlength,rgb(atom_type(i),:),style
      end do

    end if

    close(icmtx)

    return

  end subroutine write_cmtx
  
  
!===============================================================================

  subroutine write_castep_cell
! ============================

!--------------------------------------------------------------------------
!
!     This subroutine writes the configuration in cif format, suitable
!     to be read by CrystalMaker or another crystal viewer program
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use annotations
    use version_data
    use channel_numbers
    use utilities

    implicit none

    integer :: i,j,n
    integer :: ifrac1(3),ifrac2(3)
    character(len=2) :: c1(3),c2(3)
    character(len=26) :: cstring
    character(len=132) :: fname
    double precision :: symarray(3,4)
    logical :: lokay

    n = index(fileroot,'.')
    if (n>0) then
      fname = trim(wfolder)//fileroot(1:n)//'cell'
    else
      fname = trim(wfolder)//trim(fileroot)//'.cell'
    end if
    inquire(file=trim(fname),exist=lokay)
    if (lokay) then
      n = index(fname,'.',back=.true.)
      fname = 'new_'//trim(fname)
    end if
    open(icastep,file=trim(fname),form='formatted',status='unknown')

    write(icastep,'(a)') '#'
    write(icastep,'(a)') '# CASTEP cell file written by data2config, version '// &
                         trim(version_number(nversions))
    write(icastep,'(a)') '#'

    write(icastep,'(/a)') '%BLOCK LATTICE_ABC'
    write(icastep,'(3f12.6)') a,b,c
    write(icastep,'(3f12.6)') alpha,beta,gamma
    write(icastep,'(a)') '%ENDBLOCK LATTICE_ABC'

    write(icastep,'(/a)') '%BLOCK cell_constraints'
    write(icastep,'(a)') '      1   2   3'
    write(icastep,'(a)') '      0   0   0'
    write(icastep,'(a)') '%ENDBLOCK cell_constraints'

    write(icastep,'(/a)') '%BLOCK_POSITIONS_FRAC'
    do i = 1,natoms
      write(icastep,'(a2,4x,3f10.5)') element(atom_type(i)),xf(i),yf(i),zf(i)
    end do
    write(icastep,'(a)') '%ENDBLOCK POSITIONS_FRAC'

    write(icastep,'(/a)') '%BLOCK symmetry_ops'
    do  i = 1,nsym
      cstring = adjustl(coperators(i))
      call convert_coperators(c1,c2,ifrac1,ifrac2,cstring)
      call operators_to_matrix(c1,c2,ifrac1,ifrac2,symarray)
      write(icastep,'(a,i0)') '# Symmetry operator ',i
      do j = 1,4
        write(icastep,'(3f12.6)') symarray(:,j)
      end do
    end do
    write(icastep,'(a)') '%ENDBLOCK symmetry_ops'

    close(icastep)

    return

  end subroutine write_castep_cell

!===============================================================================


  subroutine write_atomeye
! ========================

!--------------------------------------------------------------------------
!
!     This subroutine writes the configuration in the atomeye cfg format
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use annotations
    use channel_numbers

    implicit none

    integer :: i,j,n,natomsuse
    character(len=132) :: text

    natomsuse = count(latomuse)

    n = index(fileroot,'.')
    if (n>0) then
      open(iatomeye,file=trim(wfolder)//fileroot(1:n)//'cfg',form='formatted',status='unknown')
    else
      open(iatomeye,file=trim(wfolder)//trim(fileroot)//'.cfg',form='formatted',status='unknown')
    endif

    write(iatomeye,'(a,i0)') 'Number of particles = ',natomsuse
    write(iatomeye,'(a)') 'A = 1.0 Angstrom'
    do i = 1,3
      do j = 1,3
        write(text,*) cell(j,i)
        text=adjustl(text) ; n = index(text,'D') ; if (n>0) text(n:n) = 'E'
        write(iatomeye,'(a,i0,a,i0,a)') 'H0(',j,',',i,') = '//trim(text)
      end do
    end do
    write(iatomeye,'(a)') '.NO_VELOCITY.'
    write(iatomeye,'(a)') 'entry_count = 3'
    
    do i = 1,natoms
      if (.not.latomuse(i)) cycle
      write(text,*) element_mass(atom_type(i))
      text=adjustl(text) ; n = index(text,'D') ; if (n>0) text(n:n) = 'E'
      write(iatomeye,'(a)') trim(text)
      write(iatomeye,'(a)') element(atom_type(i))
      write(iatomeye,*) xf(i),yf(i),zf(i)
    end do
    close(iatomeye)

    return

  end subroutine write_atomeye


!===============================================================================


  subroutine write_gasp
! =====================

!--------------------------------------------------------------------------
!
!     This subroutine writes the configuration in the format for GASP
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use annotations
    use channel_numbers

    implicit none

    integer :: i,n
    character(len=132) :: text

    n = index(fileroot,'.')
    if (n>0) then
      open(igasp,file=trim(wfolder)//fileroot(1:n)//'dat',form='formatted',status='unknown')
    else
      open(igasp,file=trim(wfolder)//trim(fileroot)//'.dat',form='formatted',status='unknown')
    endif

    if (len_trim(config_title)>0) then
      write(igasp,'(a)') trim(config_title)
    else
      write(igasp,'(a)') 'Gasp file'
    end if

    if (lwriteo) then
      do i = 1,3
        write(igasp,*) cell(:,i)
      end do    
      do i = 1,natoms
        write(igasp,'(a2,6x,3f12.6)') element(atom_type(i)),xo(i),yo(i),zo(i)
      end do
    else
      write(igasp,'(6f14.6)') a,b,c,alpha,beta,gamma
      write(igasp,'(a)') 'fractional'
      do i = 1,natoms
        write(igasp,'(a2,6x,3f12.6)') element(atom_type(i)),xf(i),yf(i),zf(i)
      end do
    end if

    close(igasp)

    return

  end subroutine write_gasp


!===============================================================================


  subroutine write_gulp
! =====================

!--------------------------------------------------------------------------
!
!     This subroutine writes the configuration in the format for GULP
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use annotations
    use version_data
    use channel_numbers

    implicit none

    integer :: i,n,nt
    character(len=132) :: buffer,buffer2

    n = index(fileroot,'.')
    if (n>0) then
      open(igulp,file=trim(wfolder)//fileroot(1:n)//'gulp',form='formatted',status='unknown')
    else
      open(igulp,file=trim(wfolder)//trim(fileroot)//'.gulp',form='formatted',status='unknown')
    end if

    buffer = '' ; buffer2 = ''
    if (allocated(norder)) then
      if (sum(norder)>0) then
        do i = 1,ntypes
          n = len_trim(buffer) + 2
          write(buffer(n:),'(a)') element(norder(i))
          n = len_trim(buffer2) + 2
          write(buffer2(n:),'(i0)') n_elements(norder(i))
        end do
      end if
    end if
    if (len_trim(buffer)==0) then
      nt = 0
      do i = -1,nelements
        if (n_elements(i)>0) then
          n = len_trim(buffer) + 2
          write(buffer(n:),'(a)') element(i)
          n = len_trim(buffer2) + 2
          write(buffer2(n:),'(i0)') n_elements(i)
          nt = nt + 1
        end if
      end do
      ntypes = nt
    end if

    write(igulp,'(a)') 'Cell'
    write(igulp,'(6f12.6)') a,b,c,alpha,beta,gamma

    if (lgulpc) then
      write(igulp,'(a,i0)') 'cartesian  ',natoms
    else
      write(igulp,'(a,i0)') 'fractional  ',natoms
    end if

    do i = 1,natoms
      if (lgulpc) then
        write(igulp,'(a,3f12.6)') element(atom_type(i)),xo(i),yo(i),zo(i)
      else
        write(igulp,'(a,3f12.6)') element(atom_type(i)),xf(i),yf(i),zf(i)
      end if
    end do

    close(igulp)
    
    return
    end subroutine write_gulp

!===============================================================================

  subroutine write_dlpoly
! =======================

!--------------------------------------------------------------------------
!
!     This subroutine writes the configuration in the DLPOLY CONFIG file
!     format
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use annotations
    use utilities
    use channel_numbers

    implicit none

    integer :: molnum,i,j,k,l,m,ierror,ibond,iangle,irigid,nbonds, &
               nangles,iline,nlines,n,nc,il,nlinks,nb,na,ntorsions,itorsion,ninversions, &
               iinversion,nrigids,iatom,jatom,katom,nrigid_bodies, &
               junk,ij,irbatom,latom,iconstraint,nconstraints,item, &
               nrigid_bodies_temp
    integer :: nlabel(4),ipad(1)
    integer, allocatable :: nrigid(:),ntemp(:),rigidtemp(:),nrigid_body(:), &
                            natoms_in_rigid_body(:),jatoms(:)
    logical :: lcontinue,lokay,ltemplate
    logical, allocatable :: lrigid(:)
    character(len=1000),allocatable :: crigid(:),ctemp(:)
    character(len=1000) :: cpad(1)
    character(len=800)  :: text,buffer,bondinfo,bondstuff,constraintinfo,constraintstuff
    character(len=4000) :: angleinfo,anglestuff,torsioninfo,torsionstuff,inversioninfo, &
                           inversionstuff
    character(len=4000),allocatable :: bond_data(:),angle_data(:),bond_info(:,:), &
                                     torsion_data(:),inversion_data(:),rigid_data(:), &
                                     constraint_data(:),constraint_info(:,:)
    character(len=4000),allocatable :: angle_info(:,:),torsion_info(:,:),inversion_info(:,:)
    character(len=4),allocatable :: rigid_pair(:)
    character(len=2) :: catom,celement(4)
    character(len=10) :: clabel(4)
    character(len=4) :: bond_pair,cb1,cb2,constraint_pair
    character(len=6) :: angle_triplet
    character(len=8) :: torsion_quartet,inversion_quartet
    double precision :: pi,rbond,dr,rmin,rmax,dx,dy,dz,dxo,dyo,dzo, &
             r,r1,r2,r3,r4,r12,r13,r23,r24,r34,tmp,tmp1,tmp2, &
             r1min,r1max,r2min,r2max,r3min,r3max,r4min,r4max,r34min,r34max, &
             dx1,dx2,dx3,dx4,dx12,dx13,dx23,dx24,dx34, &
             dy1,dy2,dy3,dy4,dy12,dy13,dy23,dy24,dy34, &
             dz1,dz2,dz3,dz4,dz12,dz13,dz23,dz24,dz34, &
             dx1o,dx2o,dx3o,dx4o,dx12o,dx13o,dx23o,dx24o,dx34o, &
             dy1o,dy2o,dy3o,dy4o,dy12o,dy13o,dy23o,dy24o,dy34o, &
             dz1o,dz2o,dz3o,dz4o,dz12o,dz13o,dz23o,dz24o,dz34o, &
             angmean,ang1,ang2,dang,angmin,angmax,ang1min,ang1max,ang2min,ang2max, &
             angle,angle1,angle2
    double precision,allocatable :: rbmin(:),rbmax(:)

    text = ''
    pi = dacos(-1.0d0)
    if (lm_title) text = adjustl(trim(config_title))
    if (lm_material) text = trim(text)//'; material = '//adjustl(trim(metadata_material))
    if (len(trim(text))==0) text = trim(text)//'; Date = '//adjustl(trim(metadata_date))

! Write the CONFIG file
    inquire(file=trim(wfolder)//'CONFIG',exist=lokay)
    if (lokay) then
      open(idlpoly,file=trim(wfolder)//'CONFIG.new',form='formatted',status='unknown')
    else
      open(idlpoly,file=trim(wfolder)//'CONFIG',form='formatted',status='new')
    end if

    write(idlpoly,'(a)') adjustl(trim(text))
    write(idlpoly,'(3i10)') 0,3,natoms
!    write(idlpoly,'(3f20.10)') cell(:,1)
!    write(idlpoly,'(3f20.10)') cell(:,2)
!    write(idlpoly,'(3f20.10)') cell(:,3)
   do i=1,3 ; write(idlpoly,'(3f20.10)') cell(i,:) ; end do

    do i = 1,natoms
    molnum = 3*(i-1)+1
    if (luselabels) then 
      if (atom_label(i)>0) then
        write(idlpoly,'(a,i0,8x,i0)') trim(element(atom_type(i))),atom_label(i),i
      else
        write(idlpoly,'(a2,8x,i0)') element(atom_type(i)),i
      end if
    else
      write(idlpoly,'(a2,8x,i0)') element(atom_type(i)),i
    end if
    write(idlpoly,'(3f20.10)') xo(i),yo(i),zo(i)
    end do
    close(idlpoly)

! Return if we don't need to construct the FIELD file
    if (.not.lfield) return

! Now we construct the FIELD file if required
    open(ifield,file=trim(fieldfile),form='formatted',status='old')
    inquire(file='template.dat',exist=ltemplate)
    if (ltemplate) open(itemplate,file='template.dat',form='formatted',status='old')

! Count number of lines in fieldfile
    ierror = 0 ; nlines = 0
    do while (ierror==0)
    read(ifield,'(a)',iostat=ierror) text
    if (ierror/=0) exit
    nlines = nlines + 1
    end do
    rewind(ifield)

! Count number of rigid terms
    ierror = 0 ; nrigids = 0
    do iline = 1,nlines
      read(ifield,'(a)',iostat=ierror) text
      text = adjustl(text)
      if (index(text,'rigid')>0) nrigids = nrigids + 1
    end do
    rewind(ifield)

! Now obtain rigid term information
    ierror = 0 ; irigid = 0
    if (nrigids>0) then
    allocate(rigid_data(nrigids))
    do iline = 1,nlines
      read(ifield,'(a)',iostat=ierror) text
      text = adjustl(text)
      if (index(text,'rigid')==0) cycle
      irigid = irigid + 1
      rigid_data(irigid) = trim(text)
    end do
    rewind(ifield)
    end if

! Count number of bonds
    ierror = 0 ; nbonds = 0
    do iline = 1,nlines
    read(ifield,'(a)',iostat=ierror) text
    text = adjustl(text)
    if (index(text,'bond')>0) nbonds = nbonds + 1
    end do
    rewind(ifield)

! Now obtain bond information
    ierror = 0 ; ibond = 0
    if (nbonds>0) then
    allocate(bond_data(nbonds),bond_info(natoms,nbonds))
    bond_info = ''
    do iline = 1,nlines
      read(ifield,'(a)',iostat=ierror) text
      text = adjustl(text)
      if (index(text,'bond')==0) cycle
      ibond = ibond + 1
      bond_data(ibond) = trim(text)
    end do
    rewind(ifield)
    end if

! Count number of constraints
    ierror = 0 ; nconstraints = 0
    do iline = 1,nlines
    read(ifield,'(a)',iostat=ierror) text
    text = adjustl(text)
    if (index(text,'constraint')>0) nconstraints = nconstraints + 1
    end do
    rewind(ifield)

! Now obtain constraint information
    ierror = 0 ; iconstraint = 0
    if (nconstraints>0) then
    allocate(constraint_data(nconstraints),constraint_info(natoms,nconstraints))
    constraint_info = ''
    do iline = 1,nlines
      read(ifield,'(a)',iostat=ierror) text
      text = adjustl(text)
      if (index(text,'constraint')==0) cycle
      iconstraint = iconstraint + 1
      constraint_data(iconstraint) = trim(text)
    end do
    rewind(ifield)
    end if

! Count number of angles
    ierror = 0 ; nangles = 0
    do iline = 1,nlines
    read(ifield,'(a)',iostat=ierror) text
    text = adjustl(text)
    if (index(text,'angle')>0) nangles = nangles + 1
    end do
    rewind(ifield)

! Now obtain angle information
    ierror = 0 ; iangle = 0
    if (nangles>0) then
    allocate(angle_data(nangles),angle_info(natoms,nangles))
    angle_info = ''
    do iline = 1,nlines
      read(ifield,'(a)',iostat=ierror) text
      text = adjustl(text)
      if (index(text,'angle')==0) cycle
      iangle = iangle + 1
      angle_data(iangle) = trim(text)
    end do
    rewind(ifield)
    end if

! Count number of torsions
    ierror = 0 ; ntorsions = 0
    do iline = 1,nlines
    read(ifield,'(a)',iostat=ierror) text
    text = adjustl(text)
    if (index(text,'torsion')>0) ntorsions = ntorsions + 1
    end do
    rewind(ifield)

! Now obtain torsion information
    ierror = 0 ; itorsion = 0
    if (ntorsions>0) then
    allocate(torsion_data(ntorsions),torsion_info(natoms,ntorsions))
    torsion_info = ''
    do iline = 1,nlines
      read(ifield,'(a)',iostat=ierror) text
      text = adjustl(text)
      if (index(text,'torsion')==0) cycle
      itorsion = itorsion + 1
      torsion_data(itorsion) = trim(text)
    end do
    rewind(ifield)
    end if

! Count number of inversions
    ierror = 0 ; ninversions = 0
    do iline = 1,nlines
    read(ifield,'(a)',iostat=ierror) text
    text = adjustl(text)
    if (index(text,'inversion')>0) ninversions = ninversions + 1
    end do
    rewind(ifield)

! Now obtain inversion information
    ierror = 0 ; iinversion = 0
    if (ninversions>0) then
    allocate(inversion_data(ninversions),inversion_info(natoms,ninversions))
    inversion_info = ''
    do iline = 1,nlines
      read(ifield,'(a)',iostat=ierror) text
      text = adjustl(text)
      if (index(text,'inversion')==0) cycle
      iinversion = iinversion + 1
      inversion_data(iinversion) = trim(text)
    end do
    rewind(ifield)
    end if

    close(ifield)

    write(6,'(i0,a)') nrigids,' rigid term(s) read from field information file'
    write(6,'(i0,a)') nbonds,' bond(s) read from field information file'
    write(6,'(i0,a)') nconstraints,' constraint(s) read from field information file'
    write(6,'(i0,a)') nangles,' angle(s) read from field information file'
    write(6,'(i0,a)') ntorsions,' torsion(s) read from field information file'
    write(6,'(i0,a)') ninversions,' inversion(s) read from field information file'

!    do irigid = 1,nrigids
!    write(6,'(a)') trim(rigid_data(irigid))
!    end do
!    do ibond = 1,nbonds
!    write(6,'(a)') trim(bond_data(ibond))
!    end do
!    do iangle = 1,nangles
!    write(6,'(a)') trim(angle_data(iangle))
!    end do
!    do itorsion = 1,ntorsions
!    write(6,'(a)') trim(torsion_data(itorsion))
!    end do
!    do iinThe order from your FIELD file
! 1,ninversions
!    write(6,'(a)') trim(inversion_data(iinversion))
!    end do

! Open the FIELDout file
   open(ifieldout,file=trim(wfolder)//'FIELDout',form='formatted',status='unknown')

!  Copy the header from the template file (called template.dat)
   if (ltemplate) then
     do
       read(itemplate,'(a)') buffer
       buffer = adjustl(buffer)
       if (buffer(1:6)=='------') exit
       write(ifieldout,'(a)') trim(buffer)
     end do
   end if

   rigid_body_if: if (nrigids>0) then


! To find rigid bodies:
! 1. Loop over all atoms, exit loop if atom is already allocated to a rigid body
! 2. Start a new rigid body, increment nrigid_bodies by one, and start a new text string
!    crigid in which the first number is the atom number. Also create an integer array
!    called natoms_in_rigid_body(nrigid_bodies) that contains the number 1 for this entry.
!    I will likely need to reallocate this array but I now know this is easy to do.
! 3. Now do a loop that ends when we add no more atoms to the rigid body. For this loop,
!    set n = natoms_in_rigid_body(nrigid_bodies) and loop over n. We can end the loop
!    when at the end n still equals natoms_in_rigid_body(nrigid_bodies).
! 4. Loop over the number of atoms within the text string. For each atom, loop over all its
!    neighbours, moving on if their logical lrigid is already true, and then loop over
!    the rigid bonds specifiers to see whether they are part of the rigid body. If they are
!    a) set lrigid to true; b) increment natoms_in_rigid_body(nrigid_bodies) by 1; c) add
!    the atom number to the end of the crigid text string.

     allocate(nrigid_body(natoms))
     allocate (lrigid(natoms))
     allocate (rigid_pair(nrigids),rbmin(nrigids),rbmax(nrigids))
     lrigid = .false.
     nrigid_body = 0
     nrigid_bodies = 0
     if (allocated(crigid)) crigid = ''
     cpad = ''
     do irigid = 1,nrigids
       buffer = adjustl(rigid_data(irigid))
       buffer = adjustl(buffer(6:))   ! After the rigid word
       bond_pair(1:2) = buffer(1:2)   ! Almost all this bit is the same as the bond part
       buffer = adjustl(buffer(3:))
       bond_pair(3:4) = buffer(1:2)
       rigid_pair(irigid) = bond_pair
       buffer = adjustl(buffer(3:))
       read(buffer,*) rbond,dr
       rbmin(irigid) = rbond - dr
       rbmax(irigid) = rbond + dr
     end do

     rigid_atom_loop_1: do iatom = 1,natoms
       if (lrigid(iatom)) cycle rigid_atom_loop_1
       nrigid_bodies = nrigid_bodies + 1   ! increment number of rigid bodies
       lrigid(iatom) = .true.
       if (allocated(nrigid)) then
         nrigid = reshape(nrigid,(/nrigid_bodies/),(/1/))
       else
         allocate(nrigid(1))
       end if
       if (allocated(natoms_in_rigid_body)) then
         natoms_in_rigid_body = reshape(natoms_in_rigid_body,(/nrigid_bodies/),(/1/))
         write(cpad(1),'(i0)') iatom
         crigid = reshape(crigid,(/nrigid_bodies/),cpad)
       else
         allocate(natoms_in_rigid_body(1))
         natoms_in_rigid_body(1) = 1
         allocate(crigid(1))
         write(crigid(1),'(i0)') iatom
       end if
       loop_over_list: do
         m = natoms_in_rigid_body(nrigid_bodies)
         if (allocated(jatoms)) deallocate(jatoms)
         allocate(jatoms(m))
         read(crigid(nrigid_bodies),*) jatoms
         loop_over_atoms_in_list: do item = 1,m
           i = jatoms(item)
           rigid_atom_loop_2: do jatom = 1,neighbours(i,0)
             j = neighbours(i,jatom)
             if (i==j) cycle rigid_atom_loop_2
             if (lrigid(j)) cycle rigid_atom_loop_2
             rigids_loop_2: do irigid = 1,nrigids
               cb1 = element(atom_type(i))//element(atom_type(j))
               cb2 = element(atom_type(j))//element(atom_type(i))
               if ((trim(cb1)/=trim(rigid_pair(irigid))).and.(trim(cb2)/=trim(rigid_pair(irigid)))) cycle rigids_loop_2
               dx = xf(i) - xf(j) + 1.5d0
               dx = dx - aint(dx) - 0.5d0
               dy = yf(i) - yf(j) + 1.5d0
               dy = dy - aint(dy) - 0.5d0
               dz = zf(i) - zf(j) + 1.5d0
               dz = dz - aint(dz) - 0.5d0
               dxo = cell(1,1)*dx + cell(1,2)*dy + cell(1,3)*dz
               dyo = cell(2,1)*dx + cell(2,2)*dy + cell(2,3)*dz
               dzo = cell(3,1)*dx + cell(3,2)*dy + cell(3,3)*dz
               r = dsqrt(dxo**2 + dyo**2 + dzo**2)
               if (r<rbmin(irigid)) cycle
               if (r>rbmax(irigid)) cycle
               natoms_in_rigid_body(nrigid_bodies) = natoms_in_rigid_body(nrigid_bodies) + 1
               n = len_trim(crigid(nrigid_bodies))+1
               write(crigid(nrigid_bodies)(n:),'(1x,i0)') j
               lrigid(j) = .true.
               lcontinue = .true.
             end do rigids_loop_2
           end do rigid_atom_loop_2
         end do loop_over_atoms_in_list
         if (m==natoms_in_rigid_body(nrigid_bodies)) exit loop_over_list
       end do loop_over_list
     end do rigid_atom_loop_1

! This bit exists because in the previous stuff we count individual atoms as
! rigid bodies, so here we need to decrement the total count
nrigid_bodies_temp = nrigid_bodies
do i = 1,nrigid_bodies
  if (natoms_in_rigid_body(i)==1) nrigid_bodies_temp = nrigid_bodies_temp - 1
end do
nrigid_bodies = nrigid_bodies_temp

  ! Write out rigid data
    write(ifieldout,'(a,1x,i0)') 'rigid',nrigid_bodies
    do i=1,nrigid_bodies
      if (natoms_in_rigid_body(i)==1) cycle
      allocate(rigidtemp(natoms_in_rigid_body(i)))
      read(crigid(i),*) rigidtemp
      write(ifieldout,'(16(i0,2x))') natoms_in_rigid_body(i),rigidtemp
      deallocate(rigidtemp)
    end do
	
  end if rigid_body_if

! Locate bonds
  do ibond = 1,nbonds
    buffer = adjustl(bond_data(ibond))
    buffer = adjustl(buffer(6:))   ! After the bond word
    n = index(buffer,' ') - 1   ! Locate the first element
    call extract_atom_name_label(buffer(1:n),celement(1),clabel(1))
    if (trim(clabel(1))=='*') then
      nlabel(1) = -2
    else if (len_trim(clabel(1))>0) then
      read(clabel(1),*,iostat=ierror) nlabel(1)
    else
      nlabel(1) = 0
    end if
    buffer = adjustl(buffer(n+1:))
    n = index(buffer,' ') - 1   ! Locate the first element
    call extract_atom_name_label(buffer(1:n),celement(2),clabel(2))
    if (trim(clabel(2))=='*') then
      nlabel(2) = -2
    else if (len_trim(clabel(2))>0) then
      read(clabel(2),*,iostat=ierror) nlabel(2)
    else
      nlabel(2) = 0
    end if
    buffer = adjustl(buffer(n+1:))
    read(buffer,*) rbond,dr
    rmin = rbond - dr
    rmax = rbond + dr
    do i = 1,natoms
      if (trim(element(atom_type(i)))/=trim(celement(1))) cycle
      if ((atom_label(i)/=nlabel(1)).and.(trim(clabel(1))/='*')) cycle
      do jatom = 1,neighbours(i,0)
        j = neighbours(i,jatom)
        if (i==j) cycle
        if (element(atom_type(j))/=celement(2)) cycle
        if ((atom_label(j)/=nlabel(2)).and.(trim(clabel(2))/='*')) cycle
        if ( (element(atom_type(j))==element(atom_type(i))) .and. (i>j) ) cycle
        dx = xf(i) - xf(j) + 1.5d0
        dx = dx - aint(dx) - 0.5d0
        dy = yf(i) - yf(j) + 1.5d0
        dy = dy - aint(dy) - 0.5d0
        dz = zf(i) - zf(j) + 1.5d0
        dz = dz - aint(dz) - 0.5d0
        dxo = cell(1,1)*dx + cell(1,2)*dy + cell(1,3)*dz
        dyo = cell(2,1)*dx + cell(2,2)*dy + cell(2,3)*dz
        dzo = cell(3,1)*dx + cell(3,2)*dy + cell(3,3)*dz
        r = dsqrt(dxo**2 + dyo**2 + dzo**2)
        if (r>rmax) cycle
        if (r<rmin) cycle
        n = len_trim(bond_info(i,ibond)) + 2
        write(bond_info(i,ibond)(n:),'(i0)') j
        n = len_trim(bond_info(i,ibond)) + 2
        bond_info(i,ibond)(n:n+1) = element(atom_type(j))
        n = len_trim(bond_info(i,ibond)) + 2
        bond_info(i,ibond)(n:n+1) = ';'
      end do
    end do
  end do

! Count the number of bonds per atom
  do i = 1,natoms
    do ibond = 1,nbonds
      n = len_trim(bond_info(i,ibond))
      nc = 0
      do j = 1,n
        if(bond_info(i,ibond)(j:j)==';') nc = nc+1
      end do
      write(bond_info(i,ibond)(n+2:),'(i0)') nc
    end do
  end do

! Count the total number of bonds
  nb = 0
  do ibond = 1,nbonds
    do i = 1,natoms
      bondinfo = adjustl(bond_info(i,ibond))
      n = index(trim(bondinfo),';',.true.)
      if (n==0) cycle
      n = n + 1
      read(bondinfo(n:),*) nlinks
      nb = nb + nlinks
    end do
  end do

! Write out the bond data
  write(ifieldout,'(a,1x,i0)') 'BONDS',nb
  do ibond = 1,nbonds
    do i = 1,natoms
      bondinfo = adjustl(bond_info(i,ibond))
      n = index(trim(bondinfo),';',.true.)
      if (n==0) cycle
      n = n + 1
      read(bondinfo(n:),*) nlinks
      n = index(bond_data(ibond),':') + 1
      bondstuff = adjustl(bond_data(ibond)(n:))
      do il = 1,nlinks
        read(bondinfo,*) j
        n = index(bondinfo,';') + 1
        bondinfo = adjustl(bondinfo(n:))
        write(ifieldout,'(a4,2x,i0,2x,i0,2x,a)') bondstuff(1:4),i,j,trim(adjustl(bondstuff(5:)))
      end do
    end do
  end do

! Locate constraints
  do iconstraint = 1,nconstraints
    buffer = adjustl(constraint_data(iconstraint))
    buffer = adjustl(buffer(12:))   ! After the constraint word
    n = index(buffer,' ') - 1   ! Locate the first element
    call extract_atom_name_label(buffer(1:n),celement(1),clabel(1))
    if (trim(clabel(1))=='*') then
      nlabel(1) = -2
    else if (len_trim(clabel(1))>0) then
      read(clabel(1),*,iostat=ierror) nlabel(1)
    else
      nlabel(1) = 0
    end if
    buffer = adjustl(buffer(n+1:))
    n = index(buffer,' ') - 1   ! Locate the first element
    call extract_atom_name_label(buffer(1:n),celement(2),clabel(2))
    if (trim(clabel(2))=='*') then
      nlabel(2) = -2
    else if (len_trim(clabel(2))>0) then
      read(clabel(2),*,iostat=ierror) nlabel(2)
    else
      nlabel(2) = 0
    end if
    buffer = adjustl(buffer(n+1:))
    read(buffer,*) rbond,dr
    rmin = rbond - dr
    rmax = rbond + dr
    do i = 1,natoms
      if (trim(element(atom_type(i)))/=trim(celement(1))) cycle
      if ((atom_label(i)/=nlabel(1)).and.(trim(clabel(1))/='*')) cycle
      do jatom = 1,neighbours(i,0)
        j = neighbours(i,jatom)
        if (i==j) cycle
        if (element(atom_type(j))/=celement(2)) cycle
        if ((atom_label(j)/=nlabel(2)).and.(trim(clabel(2))/='*')) cycle
        if ( (element(atom_type(j))==element(atom_type(i))) .and. (i>j) ) cycle
        dx = xf(i) - xf(j) + 1.5d0
        dx = dx - aint(dx) - 0.5d0
        dy = yf(i) - yf(j) + 1.5d0
        dy = dy - aint(dy) - 0.5d0
        dz = zf(i) - zf(j) + 1.5d0
        dz = dz - aint(dz) - 0.5d0
        dxo = cell(1,1)*dx + cell(1,2)*dy + cell(1,3)*dz
        dyo = cell(2,1)*dx + cell(2,2)*dy + cell(2,3)*dz
        dzo = cell(3,1)*dx + cell(3,2)*dy + cell(3,3)*dz
        r = dsqrt(dxo**2 + dyo**2 + dzo**2)
        if (r>rmax) cycle
        if (r<rmin) cycle
        n = len_trim(constraint_info(i,iconstraint)) + 2
        write(constraint_info(i,iconstraint)(n:),'(i0)') j
        n = len_trim(constraint_info(i,iconstraint)) + 2
        constraint_info(i,iconstraint)(n:n+1) = element(atom_type(j))
        n = len_trim(constraint_info(i,iconstraint)) + 2
        constraint_info(i,iconstraint)(n:n+1) = ';'
      end do
    end do
  end do

! Count the number of constraints per atom
  do i = 1,natoms
    do iconstraint = 1,nconstraints
      n = len_trim(constraint_info(i,iconstraint))
      nc = 0
      do j = 1,n
        if(constraint_info(i,iconstraint)(j:j)==';') nc = nc+1
      end do
      write(constraint_info(i,iconstraint)(n+2:),'(i0)') nc
    end do
  end do

! Count the total number of constraints
  nb = 0
  do iconstraint = 1,nconstraints
    do i = 1,natoms
      constraintinfo = adjustl(constraint_info(i,iconstraint))
      n = index(trim(constraintinfo),';',.true.)
      if (n==0) cycle
      n = n + 1
      read(constraintinfo(n:),*) nlinks
      nb = nb + nlinks
    end do
  end do

! Write out the constraint data
  write(ifieldout,'(a,1x,i0)') 'CONSTRAINTS',nb
  do iconstraint = 1,nconstraints
    do i = 1,natoms
      constraintinfo = adjustl(constraint_info(i,iconstraint))
      n = index(trim(constraintinfo),';',.true.)
      if (n==0) cycle
      n = n + 1
      read(constraintinfo(n:),*) nlinks
      n = index(constraint_data(iconstraint),':') + 1
      constraintstuff = adjustl(constraint_data(iconstraint)(n:))
      do il = 1,nlinks
        read(constraintinfo,*) j
        n = index(constraintinfo,';') + 1
        constraintinfo = adjustl(constraintinfo(n:))
        write(ifieldout,'(i0,2x,i0,2x,a)') i,j,trim(constraintstuff)
      end do
    end do
  end do

! Locate triplets
  do iangle = 1,nangles
    buffer = adjustl(angle_data(iangle))
    buffer = adjustl(buffer(7:))   ! After the angle word
    n = index(buffer,' ') - 1   ! Locate the first element
    call extract_atom_name_label(buffer(1:n),celement(1),clabel(1))
    if (trim(clabel(1))=='*') then
      nlabel(1) = -2
    else if (len_trim(clabel(1))>0) then
      read(clabel(1),*,iostat=ierror) nlabel(1)
    else
      nlabel(1) = 0
    end if
    buffer = adjustl(buffer(n+1:))
    n = index(buffer,' ') - 1   ! Locate the second element
    call extract_atom_name_label(buffer(1:n),celement(2),clabel(2))
    if (trim(clabel(2))=='*') then
      nlabel(2) = -2
    else if (len_trim(clabel(2))>0) then
      read(clabel(2),*,iostat=ierror) nlabel(2)
    else
      nlabel(2) = 0
    end if
    buffer = adjustl(buffer(n+1:))
    n = index(buffer,' ') - 1   ! Locate the third element
    call extract_atom_name_label(buffer(1:n),celement(3),clabel(3))
    if (trim(clabel(3))=='*') then
      nlabel(3) = -2
    else if (len_trim(clabel(3))>0) then
      read(clabel(3),*,iostat=ierror) nlabel(3)
    else
      nlabel(3) = 0
    end if
    buffer = adjustl(buffer(n+1:))
    read(buffer,*) r1,r2,dr,angmean,dang
    r1min = r1 - dr
    r1max = r1 + dr
    r2min = r2 - dr
    r2max = r2 + dr
    angmin = angmean - dang
    angmax = angmean + dang
    do i = 1,natoms
      if (trim(element(atom_type(i)))/=trim(celement(1))) cycle
      if ((atom_label(i)/=nlabel(1)).and.(trim(clabel(1))/='*')) cycle
      do jatom = 1,neighbours(i,0)
        j = neighbours(i,jatom)
        if (i==j) cycle
        if (trim(element(atom_type(j)))/=trim(celement(2))) cycle
        if ((atom_label(j)/=nlabel(2)).and.(trim(clabel(2))/='*')) cycle
        dx1 = xf(i) - xf(j) + 1.5d0
        dx1 = dx1 - aint(dx1) - 0.5d0
        dy1 = yf(i) - yf(j) + 1.5d0
        dy1 = dy1 - aint(dy1) - 0.5d0
        dz1 = zf(i) - zf(j) + 1.5d0
        dz1 = dz1 - aint(dz1) - 0.5d0
        dx1o = cell(1,1)*dx1 + cell(1,2)*dy1 + cell(1,3)*dz1
        dy1o = cell(2,1)*dx1 + cell(2,2)*dy1 + cell(2,3)*dz1
        dz1o = cell(3,1)*dx1 + cell(3,2)*dy1 + cell(3,3)*dz1
        r1 = dsqrt(dx1o**2 + dy1o**2 + dz1o**2)
        if (r1>r1max) cycle
        if (r1<r1min) cycle
        do katom = 1,neighbours(i,0)
          k = neighbours(i,katom)
          if ((k==i).or.(k==j)) cycle
          if ((j>k).and.(element(atom_type(j))==element(atom_type(k)))) cycle
          if (trim(element(atom_type(k)))/=trim(celement(3))) cycle
          if ((atom_label(k)/=nlabel(3)).and.(trim(clabel(3))/='*')) cycle
          dx2 = xf(i) - xf(k) + 1.5d0
          dx2 = dx2 - aint(dx2) - 0.5d0
          dy2 = yf(i) - yf(k) + 1.5d0
          dy2 = dy2 - aint(dy2) - 0.5d0
          dz2 = zf(i) - zf(k) + 1.5d0
          dz2 = dz2 - aint(dz2) - 0.5d0
          dx2o = cell(1,1)*dx2 + cell(1,2)*dy2 + cell(1,3)*dz2
          dy2o = cell(2,1)*dx2 + cell(2,2)*dy2 + cell(2,3)*dz2
          dz2o = cell(3,1)*dx2 + cell(3,2)*dy2 + cell(3,3)*dz2
          r2 = dsqrt(dx2o**2 + dy2o**2 + dz2o**2)
          if (r2>r2max) cycle
          if (r2<r2min) cycle
          dx12 = xf(j) - xf(k) + 1.5d0
          dx12 = dx12 - aint(dx12) - 0.5d0
          dy12 = yf(j) - yf(k) + 1.5d0
          dy12 = dy12 - aint(dy12) - 0.5d0
          dz12 = zf(j) - zf(k) + 1.5d0
          dz12 = dz12 - aint(dz12) - 0.5d0
          dx12o = cell(1,1)*dx12 + cell(1,2)*dy12 + cell(1,3)*dz12
          dy12o = cell(2,1)*dx12 + cell(2,2)*dy12 + cell(2,3)*dz12
          dz12o = cell(3,1)*dx12 + cell(3,2)*dy12 + cell(3,3)*dz12
          r12 = dsqrt(dx12o**2 + dy12o**2 + dz12o**2)
          tmp = (r1**2 + r2**2 - r12**2)/(2.0d0*r1*r2)
          if (tmp>=1.0d0) then
            angle = 0.0d0
          else if (tmp<=-1.0d0) then
            angle = 180.0d0
          else
            angle = 180.0d0*acos(tmp)/pi
          end if
          if ((angle<angmin).or.(angle>angmax)) cycle
          n = len_trim(angle_info(i,iangle)) + 2
          write(angle_info(i,iangle)(n:),'(i0)') j
          n = len_trim(angle_info(i,iangle)) + 2
          angle_info(i,iangle)(n:n+1) = element(atom_type(j))
          n = len_trim(angle_info(i,iangle)) + 2
          write(angle_info(i,iangle)(n:),'(i0)') k
          n = len_trim(angle_info(i,iangle)) + 2
          angle_info(i,iangle)(n:n+1) = element(atom_type(k))
          n = len_trim(angle_info(i,iangle)) + 2
          angle_info(i,iangle)(n:n+1) = ';'
        end do
      end do
    end do
  end do

! Count the number of angles per atom
  do i = 1,natoms
    do iangle = 1,nangles
      n = len_trim(angle_info(i,iangle))
      nc = 0
      do j = 1,n
        if(angle_info(i,iangle)(j:j)==';') nc = nc+1
      end do
      write(angle_info(i,iangle)(n+2:),'(i0)') nc
    end do
  end do

! Count the total number of angles
  na = 0
  do iangle = 1,nangles
    do i = 1,natoms
      angleinfo = adjustl(angle_info(i,iangle))
      n = index(trim(angleinfo),';',.true.)
      if (n==0) cycle
      n = n + 1
      read(angleinfo(n:),*) nlinks
      na = na + nlinks
    end do
  end do

! Write out the angle data
  write(ifieldout,'(a,1x,i0)') 'ANGLES',na
  do iangle = 1,nangles
    do i = 1,natoms
      angleinfo = adjustl(angle_info(i,iangle))
      n = index(trim(angleinfo),';',.true.)
      if (n==0) cycle
      n = n + 1
      read(angleinfo(n:),*) nlinks
      n = index(angle_data(iangle),':') + 1
      anglestuff = adjustl(angle_data(iangle)(n:))
      do il = 1,nlinks
        read(angleinfo,*) j,catom,k
        n = index(angleinfo,';') + 1
        angleinfo = adjustl(angleinfo(n:))
        write(ifieldout,'(a4,2x,i0,2x,i0,2x,i0,2x,a)') &
                     anglestuff(1:4),j,i,k,trim(adjustl(anglestuff(5:)))
      end do
    end do
  end do

! Locate quartets for torsional interactions
  do itorsion = 1,ntorsions
    buffer = adjustl(torsion_data(itorsion))
    buffer = adjustl(buffer(9:))   ! After the torsion word
    torsion_quartet(1:2) = buffer(1:2)
    buffer = adjustl(buffer(3:))
    torsion_quartet(3:4) = buffer(1:2)
    buffer = adjustl(buffer(3:))
    torsion_quartet(5:6) = buffer(1:2)
    buffer = adjustl(buffer(3:))
    torsion_quartet(7:8) = buffer(1:2)
    buffer = adjustl(buffer(3:))
    read(buffer,*) r1,r2,r3,r4,dr,ang1,ang2,dang
    r1min = r1 - dr
    r1max = r1 + dr
    r2min = r2 - dr
    r2max = r2 + dr
    r3min = r3 - dr
    r3max = r3 + dr
    r4min = r4 - dr
    r4max = r4 + dr
    ang1min = ang1 - dang
    ang1max = ang1 + dang
    ang2min = ang2 - dang
    ang2max = ang2 + dang
! In what follows, we consider the torsion to involve the sequence k-i-j-l    
    do i = 1,natoms
      if (element(atom_type(i))==torsion_quartet(3:4)) then
        do jatom = 1,neighbours(i,0)
          j = neighbours(i,jatom)
          if (i==j) cycle
          if (i>j) cycle
          if (element(atom_type(j))==torsion_quartet(5:6)) then
            dx2 = xf(i) - xf(j) + 1.5d0
            dx2 = dx2 - aint(dx2) - 0.5d0
            dy2 = yf(i) - yf(j) + 1.5d0
            dy2 = dy2 - aint(dy2) - 0.5d0
            dz2 = zf(i) - zf(j) + 1.5d0
            dz2 = dz2 - aint(dz2) - 0.5d0
            dx2o = cell(1,1)*dx2 + cell(1,2)*dy2 + cell(1,3)*dz2
            dy2o = cell(2,1)*dx2 + cell(2,2)*dy2 + cell(2,3)*dz2
            dz2o = cell(3,1)*dx2 + cell(3,2)*dy2 + cell(3,3)*dz2
            r2 = dsqrt(dx2o**2 + dy2o**2 + dz2o**2)
            if (r2>r2max) cycle
            if (r2<r2min) cycle
            do katom = 1,neighbours(j,0)
              k = neighbours(j,katom)
              if ((k==i).or.(k==j)) cycle
!              if (j>k) cycle
              if (element(atom_type(k))==torsion_quartet(1:2)) then
                dx1 = xf(i) - xf(k) + 1.5d0
                dx1 = dx1 - aint(dx1) - 0.5d0
                dy1 = yf(i) - yf(k) + 1.5d0
                dy1 = dy1 - aint(dy1) - 0.5d0
                dz1 = zf(i) - zf(k) + 1.5d0
                dz1 = dz1 - aint(dz1) - 0.5d0
                dx1o = cell(1,1)*dx1 + cell(1,2)*dy1 + cell(1,3)*dz1
                dy1o = cell(2,1)*dx1 + cell(2,2)*dy1 + cell(2,3)*dz1
                dz1o = cell(3,1)*dx1 + cell(3,2)*dy1 + cell(3,3)*dz1
                r1 = dsqrt(dx1o**2 + dy1o**2 + dz1o**2)
                if (r1>r1max) cycle
                if (r1<r1min) cycle
                dx13 = xf(j) - xf(k) + 1.5d0
                dx13 = dx13 - aint(dx13) - 0.5d0
                dy13 = yf(j) - yf(k) + 1.5d0
                dy13 = dy13 - aint(dy13) - 0.5d0
                dz13 = zf(j) - zf(k) + 1.5d0
                dz13 = dz13 - aint(dz13) - 0.5d0
                dx13o = cell(1,1)*dx13 + cell(1,2)*dy13 + cell(1,3)*dz13
                dy13o = cell(2,1)*dx13 + cell(2,2)*dy13 + cell(2,3)*dz13
                dz13o = cell(3,1)*dx13 + cell(3,2)*dy13 + cell(3,3)*dz13
                r13 = dsqrt(dx13o**2 + dy13o**2 + dz13o**2)
                tmp1 = (r1**2 + r2**2 - r13**2)/(2.0d0*r1*r2)
                if (tmp1>=1.0d0) then
                  angle1 = 0.0d0
                else if (tmp1<=-1.0d0) then
                  angle1 = 180.0d0
                else
                  angle1 = 180.0d0*acos(tmp)/pi
                end if
                if ((angle1<ang1min).or.(angle1>ang1max)) cycle
                !Now search for other pair
                do latom = 1,neighbours(j,0)   ! j was k
                  l = neighbours(j,latom)      ! j was k
                  if ((l==i).or.(l==j).or.(l==k)) cycle
!                  if (j>l) cycle
                  if (element(atom_type(l))==torsion_quartet(7:8)) then
                    dx3 = xf(j) - xf(l) + 1.5d0
                    dx3 = dx3 - aint(dx3) - 0.5d0
                    dy3 = yf(j) - yf(l) + 1.5d0
                    dy3 = dy3 - aint(dy3) - 0.5d0
                    dz3 = zf(j) - zf(l) + 1.5d0
                    dz3 = dz3 - aint(dz3) - 0.5d0
                    dx3o = cell(1,1)*dx3 + cell(1,2)*dy3 + cell(1,3)*dz3
                    dy3o = cell(2,1)*dx3 + cell(2,2)*dy3 + cell(2,3)*dz3
                    dz3o = cell(3,1)*dx3 + cell(3,2)*dy3 + cell(3,3)*dz3
                    r3 = dsqrt(dx3o**2 + dy3o**2 + dz3o**2)
                    if (r3>r3max) cycle
                    if (r3<r3min) cycle
                    dx24 = xf(i) - xf(l) + 1.5d0
                    dx24 = dx24 - aint(dx24) - 0.5d0
                    dy24 = yf(i) - yf(l) + 1.5d0
                    dy24 = dy24 - aint(dy24) - 0.5d0
                    dz24 = zf(i) - zf(l) + 1.5d0
                    dz24 = dz24 - aint(dz24) - 0.5d0
                    dx24o = cell(1,1)*dx24 + cell(1,2)*dy24 + cell(1,3)*dz24
                    dy24o = cell(2,1)*dx24 + cell(2,2)*dy24 + cell(2,3)*dz24
                    dz24o = cell(3,1)*dx24 + cell(3,2)*dy24 + cell(3,3)*dz24
                    r24 = dsqrt(dx24o**2 + dy24o**2 + dz24o**2)
!                    tmp2 = (r1**2 + r2**2 - r24**2)/(2.0d0*r1*r2)  ! I think this is a copy/paste error
                    tmp2 = (r2**2 + r3**2 - r24**2)/(2.0d0*r2*r3)
                    if (tmp2>=1.0d0) then
                      angle2 = 0.0d0
                    else if (tmp2<=-1.0d0) then
                      angle2 = 180.0d0
                    else
                      angle2 = 180.0d0*acos(tmp)/pi
                    end if
                    if ((angle2<ang2min).or.(angle2>ang2max)) cycle
                    !Now check the r4 distance
                    dx4 = xf(k) - xf(l) + 1.5d0
                    dx4 = dx4 - aint(dx4) - 0.5d0
                    dy4 = yf(k) - yf(l) + 1.5d0
                    dy4 = dy4 - aint(dy4) - 0.5d0
                    dz4 = zf(k) - zf(l) + 1.5d0
                    dz4 = dz4 - aint(dz4) - 0.5d0
                    dx4o = cell(1,1)*dx4 + cell(1,2)*dy4 + cell(1,3)*dz4
                    dy4o = cell(2,1)*dx4 + cell(2,2)*dy4 + cell(2,3)*dz4
                    dz4o = cell(3,1)*dx4 + cell(3,2)*dy4 + cell(3,3)*dz4
                    r4 = dsqrt(dx4o**2 + dy4o**2 + dz4o**2)
                    if (r4>r4max) cycle
                    if (r4<r4min) cycle
                    n = len_trim(torsion_info(i,itorsion)) + 2
                    write(torsion_info(i,itorsion)(n:),'(i0)') j
                    n = len_trim(torsion_info(i,itorsion)) + 2
                    torsion_info(i,itorsion)(n:n+1) = element(atom_type(j))
                    n = len_trim(torsion_info(i,itorsion)) + 2
                    write(torsion_info(i,itorsion)(n:),'(i0)') k
                    n = len_trim(torsion_info(i,itorsion)) + 2
                    torsion_info(i,itorsion)(n:n+1) = element(atom_type(k))
                    n = len_trim(torsion_info(i,itorsion)) + 2
                    write(torsion_info(i,itorsion)(n:),'(i0)') l
                    n = len_trim(torsion_info(i,itorsion)) + 2
                    torsion_info(i,itorsion)(n:n+1) = element(atom_type(l))
                    n = len_trim(torsion_info(i,itorsion)) + 2
                    torsion_info(i,itorsion)(n:n+1) = ';'
                  end if
                end do
              end if
            end do
          end if
        end do
      end if
    end do
  end do

! Count the number of torsions per atom
  do i = 1,natoms
    do itorsion = 1,ntorsions
      n = len_trim(torsion_info(i,itorsion))
      nc = 0
      do j = 1,n
        if(torsion_info(i,itorsion)(j:j)==';') nc = nc+1
      end do
      write(torsion_info(i,itorsion)(n+2:),'(i0)') nc
    end do
  end do

! Count the total number of torsions
  na = 0
  do itorsion = 1,ntorsions
    do i = 1,natoms
      torsioninfo = adjustl(torsion_info(i,itorsion))
      n = index(trim(torsioninfo),';',.true.)
      if (n==0) cycle
      n = n + 1
      read(torsioninfo(n:),*) nlinks
      na = na + nlinks
    end do
  end do

! Write out the torsion data
  write(ifieldout,'(a,1x,i0)') 'DIHEDRALS',na
  do itorsion = 1,ntorsions
    do i = 1,natoms
      torsioninfo = adjustl(torsion_info(i,itorsion))
      n = index(trim(torsioninfo),';',.true.)
      if (n==0) cycle
      n = n + 1
      read(torsioninfo(n:),*) nlinks
      n = index(torsion_data(itorsion),':') + 1
      torsionstuff = adjustl(torsion_data(itorsion)(n:))
      do il = 1,nlinks
        read(torsioninfo,*) j,catom,k,catom,l
        n = index(torsioninfo,';') + 1
        torsioninfo = adjustl(torsioninfo(n:))
        write(ifieldout,'(a4,2x,i0,2x,i0,2x,i0,2x,i0,2x,a)') &
        torsionstuff(1:4),k,i,j,l,trim(adjustl(torsionstuff(5:)))
      end do
    end do
  end do

! Locate inversions
  do iinversion = 1,ninversions
    buffer = adjustl(inversion_data(iinversion))
    buffer = adjustl(buffer(11:))   ! After the inversion word
    inversion_quartet(1:2) = buffer(1:2)
    buffer = adjustl(buffer(3:))
    inversion_quartet(3:4) = buffer(1:2)
    buffer = adjustl(buffer(3:))
    inversion_quartet(5:6) = buffer(1:2)
    buffer = adjustl(buffer(3:))
    inversion_quartet(7:8) = buffer(1:2)
    buffer = adjustl(buffer(3:))
    read(buffer,*) r2,r3,r4,dr,ang1,ang2,dang
    r2min = r2 - dr
    r2max = r2 + dr
    r3min = r3 - dr
    r3max = r3 + dr
    r4min = r4 - dr
    r4max = r4 + dr
    ang1min = ang1 - dang
    ang1max = ang1 + dang
    ang2min = ang2 - dang
    ang2max = ang2 + dang
    do i = 1,natoms
      if (element(atom_type(i))==inversion_quartet(1:2)) then
        do jatom = 1,neighbours(i,0)
          j = neighbours(i,jatom)
          if (i==j) cycle
          if (element(atom_type(j))==inversion_quartet(3:4)) then
            dx2 = xf(i) - xf(j) + 1.5d0
            dx2 = dx2 - aint(dx2) - 0.5d0
            dy2 = yf(i) - yf(j) + 1.5d0
            dy2 = dy2 - aint(dy2) - 0.5d0
            dz2 = zf(i) - zf(j) + 1.5d0
            dz2 = dz2 - aint(dz2) - 0.5d0
            dx2o = cell(1,1)*dx2 + cell(1,2)*dy2 + cell(1,3)*dz2
            dy2o = cell(2,1)*dx2 + cell(2,2)*dy2 + cell(2,3)*dz2
            dz2o = cell(3,1)*dx2 + cell(3,2)*dy2 + cell(3,3)*dz2
            r2 = dsqrt(dx2o**2 + dy2o**2 + dz2o**2)
            if (r2>r2max) cycle
            if (r2<r2min) cycle
            do katom = 1,neighbours(i,0)
              k = neighbours(i,katom)
              if ((k==i).or.(k==j)) cycle
              if ((atom_type(j)==atom_type(k)).and.(j>k)) cycle
              if (element(atom_type(k))==inversion_quartet(5:6)) then
                dx3 = xf(i) - xf(k) + 1.5d0
                dx3 = dx3 - aint(dx3) - 0.5d0
                dy3 = yf(i) - yf(k) + 1.5d0
                dy3 = dy3 - aint(dy3) - 0.5d0
                dz3 = zf(i) - zf(k) + 1.5d0
                dz3 = dz3 - aint(dz3) - 0.5d0
                dx3o = cell(1,1)*dx3 + cell(1,2)*dy3 + cell(1,3)*dz3
                dy3o = cell(2,1)*dx3 + cell(2,2)*dy3 + cell(2,3)*dz3
                dz3o = cell(3,1)*dx3 + cell(3,2)*dy3 + cell(3,3)*dz3
                r3 = dsqrt(dx3o**2 + dy3o**2 + dz3o**2)
                if (r3>r3max) cycle
                if (r3<r3min) cycle
                dx23 = xf(j) - xf(k) + 1.5d0
                dx23 = dx23 - aint(dx23) - 0.5d0
                dy23 = yf(j) - yf(k) + 1.5d0
                dy23 = dy23 - aint(dy23) - 0.5d0
                dz23 = zf(j) - zf(k) + 1.5d0
                dz23 = dz23 - aint(dz23) - 0.5d0
                dx23o = cell(1,1)*dx23 + cell(1,2)*dy23 + cell(1,3)*dz23
                dy23o = cell(2,1)*dx23 + cell(2,2)*dy23 + cell(2,3)*dz23
                dz23o = cell(3,1)*dx23 + cell(3,2)*dy23 + cell(3,3)*dz23
                r23 = dsqrt(dx23o**2 + dy23o**2 + dz23o**2)
                tmp1 = (r2**2 + r3**2 - r23**2)/(2.0d0*r2*r3)
                if (tmp1>=1.0d0) then
                  angle1 = 0.0d0
                else if (tmp1<=-1.0d0) then
                  angle1 = 180.0d0
                else
                  angle1 = 180.0d0*acos(tmp1)/pi
                end if
                if ((angle1<ang1min).or.(angle1>ang1max)) cycle
                !Now search for other pair
                do latom = 1,neighbours(i,0)
                  l = neighbours(i,latom)
                  if ((l==i).or.(l==j).or.(l==k)) cycle
                  if ((atom_type(j)==atom_type(l)).and.(j>l)) cycle
                  if ((atom_type(k)==atom_type(l)).and.(k>l)) cycle
                  if (element(atom_type(l))==inversion_quartet(7:8)) then
                    dx4 = xf(i) - xf(l) + 1.5d0
                    dx4 = dx4 - aint(dx4) - 0.5d0
                    dy4 = yf(i) - yf(l) + 1.5d0
                    dy4 = dy4 - aint(dy4) - 0.5d0
                    dz4 = zf(i) - zf(l) + 1.5d0
                    dz4 = dz4 - aint(dz4) - 0.5d0
                    dx4o = cell(1,1)*dx4 + cell(1,2)*dy4 + cell(1,3)*dz4
                    dy4o = cell(2,1)*dx4 + cell(2,2)*dy4 + cell(2,3)*dz4
                    dz4o = cell(3,1)*dx4 + cell(3,2)*dy4 + cell(3,3)*dz4
                    r4 = dsqrt(dx4o**2 + dy4o**2 + dz4o**2)
                    if (r4>r4max) cycle
                    if (r4<r4min) cycle
                    dx24 = xf(j) - xf(l) + 1.5d0
                    dx24 = dx24 - aint(dx24) - 0.5d0
                    dy24 = yf(j) - yf(l) + 1.5d0
                    dy24 = dy24 - aint(dy24) - 0.5d0
                    dz24 = zf(j) - zf(l) + 1.5d0
                    dz24 = dz24 - aint(dz24) - 0.5d0
                    dx24o = cell(1,1)*dx24 + cell(1,2)*dy24 + cell(1,3)*dz24
                    dy24o = cell(2,1)*dx24 + cell(2,2)*dy24 + cell(2,3)*dz24
                    dz24o = cell(3,1)*dx24 + cell(3,2)*dy24 + cell(3,3)*dz24
                    r24 = dsqrt(dx24o**2 + dy24o**2 + dz24o**2)
                    tmp2 = (r2**2 + r4**2 - r24**2)/(2.0d0*r2*r4)
                    if (tmp2>=1.0d0) then
                      angle2 = 0.0d0
                    else if (tmp2<=-1.0d0) then
                      angle2 = 180.0d0
                    else
                      angle2 = 180.0d0*acos(tmp2)/pi
                    end if
                    if ((angle2<ang2min).or.(angle2>ang2max)) cycle
                    !Now check the r34 distance
                    dx34 = xf(k) - xf(l) + 1.5d0
                    dx34 = dx34 - aint(dx34) - 0.5d0
                    dy34 = yf(k) - yf(l) + 1.5d0
                    dy34 = dy34 - aint(dy34) - 0.5d0
                    dz34 = zf(k) - zf(l) + 1.5d0
                    dz34 = dz34 - aint(dz34) - 0.5d0
                    dx34o = cell(1,1)*dx34 + cell(1,2)*dy34 + cell(1,3)*dz34
                    dy34o = cell(2,1)*dx34 + cell(2,2)*dy34 + cell(2,3)*dz34
                    dz34o = cell(3,1)*dx34 + cell(3,2)*dy34 + cell(3,3)*dz34
                    r34 = dsqrt(dx34o**2 + dy34o**2 + dz34o**2)
                    r34min = r34 - dr
                    r34max = r34 + dr
                    if (r34>r34max) cycle
                    if (r34<r34min) cycle
                    n = len_trim(inversion_info(i,iinversion)) + 2
                    write(inversion_info(i,iinversion)(n:),'(i0)') j
                    n = len_trim(inversion_info(i,iinversion)) + 2
                    inversion_info(i,iinversion)(n:n+1) = element(atom_type(j))
                    n = len_trim(inversion_info(i,iinversion)) + 2
                    write(inversion_info(i,iinversion)(n:),'(i0)') k
                    n = len_trim(inversion_info(i,iinversion)) + 2
                    inversion_info(i,iinversion)(n:n+1) = element(atom_type(k))
                    n = len_trim(inversion_info(i,iinversion)) + 2
                    write(inversion_info(i,iinversion)(n:),'(i0)') l
                    n = len_trim(inversion_info(i,iinversion)) + 2
                    inversion_info(i,iinversion)(n:n+1) = element(atom_type(l))
                    n = len_trim(inversion_info(i,iinversion)) + 2
                    inversion_info(i,iinversion)(n:n+1) = ';'
                  end if
                end do
              end if
            end do
          end if
        end do
      end if
    end do
  end do

! Count the number of inversions per atom
  do i = 1,natoms
    do iinversion = 1,ninversions
      n = len_trim(inversion_info(i,iinversion))
      nc = 0
      do j = 1,n
        if(inversion_info(i,iinversion)(j:j)==';') nc = nc+1
      end do
      write(inversion_info(i,iinversion)(n+2:),'(i0)') nc
    end do
  end do

! Count the total number of inversions
  na = 0
  do iinversion = 1,ninversions
    do i = 1,natoms
      inversioninfo = adjustl(inversion_info(i,iinversion))
      n = index(trim(inversioninfo),';',.true.)
      if (n==0) cycle
      n = n + 1
      read(inversioninfo(n:),*) nlinks
      na = na + nlinks
    end do
  end do

! Write out the inversion data
  write(ifieldout,'(a,1x,i0)') 'INVERSIONS',na
  do iinversion = 1,ninversions
    do i = 1,natoms
      inversioninfo = adjustl(inversion_info(i,iinversion))
      n = index(trim(inversioninfo),';',.true.)
      if (n==0) cycle
      n = n + 1
      read(inversioninfo(n:),*) nlinks
      n = index(inversion_data(iinversion),':') + 1
      inversionstuff = adjustl(inversion_data(iinversion)(n:))
      do il = 1,nlinks
        read(inversioninfo,*) j,catom,k,catom,l
        n = index(inversioninfo,';') + 1
        inversioninfo = adjustl(inversioninfo(n:))
        write(ifieldout,'(a4,2x,i0,2x,i0,2x,i0,2x,i0,2x,a)') &
        inversionstuff(1:4),i,j,k,l,trim(adjustl(inversionstuff(5:)))
      end do
    end do
  end do

  if (ltemplate) then
     do
       read(itemplate,'(a)') buffer
       buffer = adjustl(buffer)
       write(ifieldout,'(a)') trim(buffer)
       if (buffer(1:5)=='CLOSE') exit
     end do
  end if

  close(ifieldout)

  if (allocated(crigid)) deallocate(crigid)
!    if (allocated(rigid_molecule)) deallocate(rigid_molecule)
!    if (allocated(rigid_molecule_out)) deallocate(rigid_molecule_out)
  if (allocated(bond_data)) deallocate(bond_data)
  if (allocated(constraint_data)) deallocate(constraint_data)
  if (allocated(angle_data)) deallocate(angle_data)
  if (allocated(bond_info)) deallocate(bond_info)
  if (allocated(constraint_info)) deallocate(constraint_info)
  if (allocated(torsion_data)) deallocate(torsion_data)
  if (allocated(inversion_data)) deallocate(inversion_data)
  if (allocated(rigid_data)) deallocate(rigid_data)
  if (allocated(angle_info)) deallocate(angle_info)
  if (allocated(torsion_info)) deallocate(torsion_info)
  if (allocated(inversion_info)) deallocate(inversion_info)

  return

  end subroutine write_dlpoly

!===============================================================================

  subroutine write_rmc3
! =====================

!--------------------------------------------------------------------------
!
!     This subroutine writes the configuration in the RMC v3 format
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use annotations
    use channel_numbers

    implicit none

    integer :: i,j,n,ielc
    double precision :: hcell(3,3)

    n = index(fileroot,'.')
    if (n>0) then
      open(irmc3,file=trim(wfolder)//fileroot(1:n)//'cfg',form='formatted',status='unknown')
    else
      open(irmc3,file=trim(wfolder)//trim(fileroot)//'.cfg',form='formatted',status='unknown')
    end if

    do i = 1,3
    do j = 1,3
      hcell(i,j) = cell(i,j)/2.0d0
    end do
    end do

    if (.not.lm_title) then
     write(6,'(a,$)') 'Please give title for configuration: '
     read(5,'(a)') config_title
    end if

    write(irmc3,'(a)') ' (Version 3 format configuration file)'
!    write(irmc3,'(a)') ' Crystal configuration'
    write(irmc3,'(a)') trim(config_title)
    write(irmc3,*)
    write(irmc3,'(a)') '          0         0         0 moves generated, tried, accepted'
    write(irmc3,'(a)') '          0                     configurations saved'
    write(irmc3,*)
    ntypes = 0
    do i = -1,nelements
    if (n_elements(i)>0) ntypes = ntypes + 1
    end do
    write(irmc3,'(i11,1x,a)') natoms,'molecules (of all types)'
    write(irmc3,'(i11,1x,a)') ntypes,'types of molecule'
    write(irmc3,'(10x,a)') '1 is the largest number of atoms in a molecule'
    write(irmc3,'(10x,a)') '0 Euler angles are provided'
    write(irmc3,*)
    write(irmc3,'(10x,a)') 'F (Box is not truncated octahedral)'
    write(irmc3,'(12x,a)') 'Defining vectors are:'
    write(irmc3,'(11x,3f11.6)') hcell(:,1)
    write(irmc3,'(11x,3f11.6)') hcell(:,2)
    write(irmc3,'(11x,3f11.6)') hcell(:,3)
    write(irmc3,*)

    write(6,'(a)') 'Please note that atoms in v3 configuration file are given in order:'
write(6,*) norder
write(6,*) numoftype 

    ielc = -1
    do i = 1,ntypes
      if (lorder) then
        ielc = norder(i)
      else
        do while (n_elements(ielc)==0)
          ielc = ielc + 1
        end do
      end if
      write(irmc3,'(i11,1x,a,i3)') n_elements(ielc),'molecules of type',i
      write(irmc3,'(10x,a)') '1 atomic sites'
      write(irmc3,'(13x,a)') '0.000000   0.000000   0.000000'
      write(irmc3,*)
      write(6,'(i0,2a)') n_elements(ielc),' atoms of element ',trim(element(ielc))
      ielc = ielc + 1
    end do
    do i = 1,natoms
      write(irmc3,'(3f11.7)') 2.0d0*xf(i)-1.0d0, 2.0d0*yf(i)-1.0d0, 2.0d0*zf(i)-1.0d0
    end do

    return

  end subroutine write_rmc3

!===============================================================================

  subroutine write_rmc6f
! ======================

!--------------------------------------------------------------------------
!
!     This subroutine writes the configuration in the RMC v6f format
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use annotations
    use version_data
    use channel_numbers

    implicit none

    integer :: i,n,nt
    character(len=132) :: buffer,buffer2,fname
    logical :: lokay

    n = index(fileroot,'.')
    if (n>0) then
      fname = trim(wfolder)//fileroot(1:n)//'rmc6f'
    else
      fname = trim(wfolder)//trim(fileroot)//'.rmc6f'
    end if
    inquire(file=trim(fname),exist=lokay)
    if (lokay) then
      n = index(fname,'.',back=.true.)
      fname = fname(1:n-1)//'_new'//trim(fname(n:))
    end if
    open(irmc6f,file=trim(fname),form='formatted',status='unknown')

    write(irmc6f,'(a)') '(Version 6f format configuration file)'
    write(irmc6f,'(a)') '(Generated by data2config version '//trim(version_number(nversions))//')'
    if (len_trim(metadata_owner)>0) &
             write(irmc6f,'(a)') 'Metadata owner:       '// trim(metadata_owner)
    if (len_trim(metadata_affiliation)>0) &
             write(irmc6f,'(a)') 'Metadata affiliation: '// trim(metadata_affiliation)
    if (len_trim(metadata_material)>0) &
             write(irmc6f,'(a)') 'Metadata material:    '// trim(metadata_material)
    if (len_trim(metadata_phase)>0) &
             write(irmc6f,'(a)') 'Metadata phase:       '// trim(metadata_phase)
    if (len_trim(metadata_formula)>0) &
             write(irmc6f,'(a)') 'Metadata formula:     '// trim(metadata_formula)
    if (len_trim(metadata_title)>0) &
             write(irmc6f,'(a)') 'Metadata title:       '// trim(metadata_title)
    if (len_trim(metadata_purpose)>0) &
             write(irmc6f,'(a)') 'Metadata purpose:     '// trim(metadata_purpose)
    if (len_trim(metadata_keywords)>0) &
             write(irmc6f,'(a)') 'Metadata keywords:    '// trim(metadata_keywords)
    if (len_trim(metadata_temperature)>0) &
             write(irmc6f,'(a)') 'Metadata temperature: '// trim(metadata_temperature)
    if (len_trim(metadata_pressure)>0) &
             write(irmc6f,'(a)') 'Metadata pressure     '// trim(metadata_pressure)
    if (len_trim(metadata_note)>0) &
             write(irmc6f,'(a)') 'Metadata note:        '// trim(metadata_note)
    if (len_trim(metadata_comment)>0) &
             write(irmc6f,'(a)') 'Metadata comment:     '// trim(metadata_comment)
    if (len_trim(metadata_date)>0) &
             write(irmc6f,'(a)') 'Metadata date:        '// trim(metadata_date)

    buffer = '' ; buffer2 = ''
    if (allocated(norder)) then
      if (sum(norder)>0) then
        do i = 1,ntypes
          n = len_trim(buffer) + 2
          write(buffer(n:),'(a)') element(norder(i))
          n = len_trim(buffer2) + 2
          write(buffer2(n:),'(i0)') n_elements(norder(i))
        end do
      end if
    end if
    if (len_trim(buffer)==0) then
      nt = 0
      do i = -1,nelements
        if (n_elements(i)>0) then
          n = len_trim(buffer) + 2
          write(buffer(n:),'(a)') element(i)
          n = len_trim(buffer2) + 2
          write(buffer2(n:),'(i0)') n_elements(i)
          nt = nt + 1
        end if
      end do
      ntypes = nt
    end if

    write(irmc6f,'(a,i0)') 'Number of types of atoms:   ',ntypes
    write(irmc6f,'(a)')    'Atom types present:         '//adjustl(trim(buffer))
    write(irmc6f,'(a)')    'Number of each atom type:   '//adjustl(trim(buffer2))
    write(irmc6f,'(a,i0)') 'Number of moves generated:           ',moves_generated
    write(irmc6f,'(a,i0)') 'Number of moves tried:               ',moves_tried
    write(irmc6f,'(a,i0)') 'Number of moves accepted:            ',moves_accepted
    write(irmc6f,'(a,i0)') 'Number of prior configuration saves: ',prior_saves
    write(irmc6f,'(a,i0)') 'Number of atoms:                     ',natoms
    ! >>>>>>>>>>>>>> Yuanpeng -> The following change is to >>>>>>
    ! consider both the supercell information read in from the >>>
    ! input '.rmc6f' file and that from the arguments. >>>>>>>>>>>
    if (input_rmc6f) then
    	write(irmc6f,'(a,i0,1x,i0,1x,i0)') &
         'Supercell dimensions: ',ncell_rmc6f(1)*ncell(1), &
          ncell_rmc6f(2)*ncell(2), &
          ncell_rmc6f(3)*ncell(3)
    else
    	write(irmc6f,'(a,i0,1x,i0,1x,i0)') &
            'Supercell dimensions:                ',ncell
    end if
    ! <<<<<<<<<<<<<< Yuanpeng finishes here <<<<<<<<<<<<<<<<<<<<<<
    if (nspins>0) write(irmc6f,'(a,i0)') 'Number of spins:                     ',nspins
    write(irmc6f,'(a,f12.6)') 'Number density (Ang^-3):             ',density

    write(irmc6f,'(a,6f12.6)') 'Cell (Ang/deg): ',a,b,c,alpha,beta,gamma
    write(irmc6f,'(a)') 'Lattice vectors (Ang):'
    do i = 1,3
    write(irmc6f,'(3f12.6)') cell(:,i)
    end do
    write(irmc6f,'(a)') 'Atoms:'
    do i = 1,natoms
    if (lmag.and.(i<=nspins)) then
      write(irmc6f,'(i0,1x,a,1x,3f12.6,1x,i0,3(1x,i0),3x,a,2x,3f12.6)') i,element(atom_type(i)), &
         xf(i),yf(i),zf(i),reference_number(i),reference_cell(i,:),'M:',spin(i,:)
    else
      if (atom_label(i)>0) then
        write(irmc6f,'(i0,1x,a,1x,a,i0,a,1x,3f12.6,1x,i0,3(1x,i0))') i,element(atom_type(i)), &
           '[',atom_label(i),']',xf(i),yf(i),zf(i),reference_number(i),reference_cell(i,:)
      else
        write(irmc6f,'(i0,1x,a,1x,3f12.6,1x,i0,3(1x,i0))') i,element(atom_type(i)), &
           xf(i),yf(i),zf(i),reference_number(i),reference_cell(i,:)
      end if
    end if
    end do

    close(irmc6f)
    
!    if (lpotentials) call write_rmc_potentials_file

    return
    end subroutine write_rmc6f

!===============================================================================

  subroutine write_rmc6o
! ======================

!--------------------------------------------------------------------------
!
!     This subroutine writes the configuration in the RMC v6o format
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use annotations
    use version_data
    use channel_numbers

    implicit none

    integer :: i,n,nt
    character(len=132) :: buffer,buffer2

    n = index(fileroot,'.')
    if (n>0) then
      open(irmc6o,file=trim(wfolder)//fileroot(1:n)//'rmc6o',form='formatted',status='unknown')
    else
      open(irmc6o,file=trim(wfolder)//trim(fileroot)//'.rmc6o',form='formatted',status='unknown')
    end if

    write(irmc6o,'(a)') '(Version 6o format configuration file)'
    write(irmc6o,'(a)') '(Generated by data2config version '//trim(version_number(nversions))//')'
    if (len_trim(metadata_owner)>0) &
             write(irmc6o,'(a)') 'Metadata owner:       '// trim(metadata_owner)
    if (len_trim(metadata_affiliation)>0) &
             write(irmc6o,'(a)') 'Metadata affiliation: '// trim(metadata_affiliation)
    if (len_trim(metadata_material)>0) &
             write(irmc6o,'(a)') 'Metadata material:    '// trim(metadata_material)
    if (len_trim(metadata_phase)>0) &
             write(irmc6o,'(a)') 'Metadata phase:       '// trim(metadata_phase)
    if (len_trim(metadata_formula)>0) &
             write(irmc6o,'(a)') 'Metadata formula:     '// trim(metadata_formula)
    if (len_trim(metadata_title)>0) &
             write(irmc6o,'(a)') 'Metadata title:       '// trim(metadata_title)
    if (len_trim(metadata_purpose)>0) &
             write(irmc6o,'(a)') 'Metadata purpose:     '// trim(metadata_purpose)
    if (len_trim(metadata_keywords)>0) &
             write(irmc6o,'(a)') 'Metadata keywords:    '// trim(metadata_keywords)
    if (len_trim(metadata_temperature)>0) &
             write(irmc6o,'(a)') 'Metadata temperature: '// trim(metadata_temperature)
    if (len_trim(metadata_pressure)>0) &
             write(irmc6o,'(a)') 'Metadata pressure     '// trim(metadata_pressure)
    if (len_trim(metadata_note)>0) &
             write(irmc6o,'(a)') 'Metadata note:        '// trim(metadata_note)
    if (len_trim(metadata_comment)>0) &
             write(irmc6o,'(a)') 'Metadata comment:     '// trim(metadata_comment)
    if (len_trim(metadata_date)>0) &
             write(irmc6o,'(a)') 'Metadata date:        '// trim(metadata_date)

    buffer = '' ; buffer2 = ''
    if (allocated(norder)) then
      if (sum(norder)>0) then
        do i = 1,ntypes
          n = len_trim(buffer) + 2
          write(buffer(n:),'(a)') element(norder(i))
          n = len_trim(buffer2) + 2
          write(buffer2(n:),'(i0)') n_elements(norder(i))
        end do
      end if
    end if
    if (len_trim(buffer)==0) then
      nt = 0
      do i = -1,nelements
        if (n_elements(i)>0) then
          n = len_trim(buffer) + 2
          write(buffer(n:),'(a)') element(i)
          n = len_trim(buffer2) + 2
          write(buffer2(n:),'(i0)') n_elements(i)
          nt = nt + 1
        end if
      end do
    end if
    write(irmc6o,'(a,i0)') 'Number of types of atoms:   ',nt
    write(irmc6o,'(a)')    'Atom types present:         '//adjustl(trim(buffer))
    write(irmc6o,'(a)')    'Number of each atom type:   '//adjustl(trim(buffer2))
    write(irmc6o,'(a,i0)') 'Number of moves generated:           ',moves_generated
    write(irmc6o,'(a,i0)') 'Number of moves tried:               ',moves_tried
    write(irmc6o,'(a,i0)') 'Number of moves accepted:            ',moves_accepted
    write(irmc6o,'(a,i0)') 'Number of prior configuration saves: ',prior_saves
    write(irmc6o,'(a,i0)') 'Number of atoms:                     ',natoms
    write(irmc6o,'(a,f12.6)') 'Number density (Ang^-3):             ',density

    write(irmc6o,'(a)') 'Lattice vectors (Ang):'
    do i = 1,3
    write(irmc6o,'(3f12.6)') cell(:,i)
    end do
    write(irmc6o,'(a)') 'Atoms:'
    do i = 1,natoms
      if (atom_label(i)>0) then
        write(irmc6o,'(i0,1x,a,1x,a,i0,a,1x,3f12.6,1x,i0,3(1x,i0))') i,element(atom_type(i)), &
            '[',atom_label(i),']',xo(i),yo(i),zo(i),reference_number(i),reference_cell(i,:)
      else
        write(irmc6o,'(i0,1x,a,1x,3f12.6,1x,i0,3(1x,i0))') i,element(atom_type(i)), &
           xo(i),yo(i),zo(i),reference_number(i),reference_cell(i,:)
      end if
    end do

    close(irmc6o)

    return
    end subroutine write_rmc6o

!===============================================================================

  subroutine write_his6f
! ======================

!--------------------------------------------------------------------------
!
!     This subroutine writes the history file in the RMC v6f format
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use annotations
    use histogram_commons
    use channel_numbers

   implicit none

    integer :: n,i

    n = index(fileroot,'.')
    if (n>0) then
      open(ihis6f,file=trim(wfolder)//fileroot(1:n)//'his6f',form='formatted',status='unknown')
    else
      open(ihis6f,file=trim(wfolder)//trim(fileroot)//'.his6f',form='formatted',status='unknown')
    end if

    write(ihis6f,'(a)') 'RMCProfile v6f intermediate (histogram) file'

    if (lm_title)    write(ihis6f,'(a)') 'Metadata tile:      '//trim(config_title)
    if (lm_owner)    write(ihis6f,'(a)') 'Metadata owner:     '//trim(metadata_owner)
    write(ihis6f,'(a)') 'Metadata date:      '//trim(metadata_date)
    if (lm_material) write(ihis6f,'(a)') 'Metadata material:  '//trim(metadata_material)
    if (lm_comment)  write(ihis6f,'(a)') 'Metadata comment:   '//trim(metadata_comment)
    if (lm_source)   write(ihis6f,'(a)') 'Metadata source:    '//trim(metadata_source)

    write(ihis6f,'(a,i0)') 'Number of moves generated:           ',moves_generated
    write(ihis6f,'(a,i0)') 'Number of moves tried:               ',moves_tried
    write(ihis6f,'(a,i0)') 'Number of moves accepted:            ',moves_accepted
    write(ihis6f,'(a,i0)') 'Number of prior configuration saves: ',prior_saves
    write(ihis6f,'(a,i0)') 'Number of atoms:                     ',natoms
    write(ihis6f,'(a,f12.6)') 'Number density (Ang^-3):             ',density

    write(ihis6f,'(a,6f12.6)') 'Cell (Ang/deg): ',a,b,c,alpha,beta,gamma
    write(ihis6f,'(a)') 'Lattice vectors (Ang):'
    do i = 1,3
    write(ihis6f,'(3f12.6)') cell(:,i)
    end do
    write(ihis6f,'(a)') 'Atoms (fractional coordinates):'
    do i = 1,natoms
    write(ihis6f,'(i6,1x,a4,3f12.6)') i,element(atom_type(i)),xf(i),yf(i),zf(i)
    end do

    write(ihis6f,'(a,i0)') 'Number of points in pair distribution functions: ',nr
    write(ihis6f,*) 'Step size in pair distribution function: ',dr
    write(ihis6f,*)
    write(ihis6f,*) 'step    ',(apairs(i)//'    ',i=1,npar)
    do i = 1,nr
    write(ihis6f,*) i,histogram(i,:)
    end do
    if (ncoord_0>0) then
    write(ihis6f,*)
    write(ihis6f,*) 'typec_0:'
    write(ihis6f,*) typec_0
    write(ihis6f,*) 'typen_0:'
    write(ihis6f,*) typen_0
    write(ihis6f,*) 'rcoord_0(1):'
    write(ihis6f,*) rcoord_0(1,:)
    write(ihis6f,*) 'rcoord_0(2):'
    write(ihis6f,*) rcoord_0(2,:)
    write(ihis6f,*) 'nneigh:'
    do i = 1,natoms
      write(ihis6f,*) nneigh(i,:)
    end do
    end if

    return

  end subroutine write_his6f

!===============================================================================

  subroutine write_his6o
! ======================

!--------------------------------------------------------------------------
!
!     This subroutine writes the history file in the RMC v6o format
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use annotations
    use histogram_commons
    use channel_numbers

    implicit none

    integer :: n,i

    n = index(fileroot,'.')
    if (n>0) then
      open(ihis6o,file=trim(wfolder)//fileroot(1:n)//'his6o',form='formatted',status='unknown')
    else
      open(ihis6o,file=trim(wfolder)//trim(fileroot)//'.his6o',form='formatted',status='unknown')
    end if

    write(6,'(a)') 'RMCProfile v6o intermediate (histogram) file'

    if (lm_title)    write(ihis6o,'(a)') 'Metadata tile:      '//trim(config_title)
    if (lm_owner)    write(ihis6o,'(a)') 'Metadata owner:     '//trim(metadata_owner)
    write(ihis6o,'(a)') 'Metadata date:      '//trim(metadata_date)
    if (lm_material) write(ihis6o,'(a)') 'Metadata material:  '//trim(metadata_material)
    if (lm_comment)  write(ihis6o,'(a)') 'Metadata comment:   '//trim(metadata_comment)
    if (lm_source)   write(ihis6o,'(a)') 'Metadata source:    '//trim(metadata_source)

    write(ihis6o,'(a,i0)') 'Number of moves generated:           ',moves_generated
    write(ihis6o,'(a,i0)') 'Number of moves tried:               ',moves_tried
    write(ihis6o,'(a,i0)') 'Number of moves accepted:            ',moves_accepted
    write(ihis6o,'(a,i0)') 'Number of prior configuration saves: ',prior_saves
    write(ihis6o,'(a,i0)') 'Number of atoms:                     ',natoms
    write(ihis6o,'(a,f12.6)') 'Number density (Ang^-3):             ',density

    write(ihis6o,'(a,6f12.6)') 'Cell (Ang/deg): ',a,b,c,alpha,beta,gamma
    write(ihis6o,'(a)') 'Lattice vectors (Ang):'
    do i = 1,3
    write(ihis6o,'(3f12.6)') cell(:,i)
    end do
    write(ihis6o,'(a)') 'Atoms (cartesian coordinates/Ang):'
    do i = 1,natoms
    write(ihis6o,'(i6,1x,a4,3f12.6)') i,element(atom_type(i)),xo(i),yo(i),zo(i)
    end do

    write(ihis6o,'(a,i0)') 'Number of points in pair distribution functions: ',nr
    write(ihis6o,*) 'Step size in pair distribution function: ',dr
    write(ihis6o,*)
    write(ihis6o,*) 'step    ',(apairs(i)//'    ',i=1,npar)
    do i = 1,nr
    write(ihis6o,*) i,histogram(i,:)
    end do
    if (ncoord_0>0) then
    write(ihis6o,*)
    write(ihis6o,*) 'typec_0:'
    write(ihis6o,*) typec_0
    write(ihis6o,*) 'typen_0:'
    write(ihis6o,*) typen_0
    write(ihis6o,*) 'rcoord_0(1):'
    write(ihis6o,*) rcoord_0(1,:)
    write(ihis6o,*) 'rcoord_0(2):'
    write(ihis6o,*) rcoord_0(2,:)
    write(ihis6o,*) 'nneigh:'
    do i = 1,natoms
      write(ihis6o,*) nneigh(i,:)
    end do
    end if

    close(ihis6o)
    return

  end subroutine write_his6o

!===============================================================================

  subroutine write_crystal
! ========================

    use structure_data
    use arguments
    use annotations
    use channel_numbers

    implicit none

    integer :: i,j,n,jatoms,iel

    if (.not.lm_title) then
     write(6,'(a)', advance="no") 'Please give title for configuration: '
     read(5,'(a)') config_title
    end if

! Sort numbers of atoms of each type
    if (.not.lorder) then
      ntypes = 0
      do iel = -1,nelements
        if (n_elements(iel)>0) then
          ntypes = ntypes + 1
        end if
      end do
      if (.not.allocated(numoftype)) allocate(numoftype(ntypes))
      ntypes = 0
      do iel = -1,nelements
        if (n_elements(iel)>0) then
          ntypes = ntypes + 1
          numoftype(ntypes) = n_elements(iel)
          write(6,'(i0,2a)') n_elements(iel),' atoms of element ',trim(element(iel))
        end if
      end do
    end if

    n = index(fileroot,'.')
    if (n>0) then
      open(icrystal,file=trim(wfolder)//fileroot(1:n)//'cfgcom',form='formatted',status='unknown')
    else
      open(icrystal,file=trim(wfolder)//trim(fileroot)//'.cfgcom',form='formatted',status='unknown')
    end if

    write(icrystal,*) 1,1,1
    do i = 1,3
      write(icrystal,'(3f12.6)') cell(:,i)
    end do
    write(icrystal,*) ntypes
      jatoms = 0
      do i = 1,ntypes
        n = numoftype(i)
        write(icrystal,*) n
        do j = 1,n
          jatoms = jatoms + 1
          write(icrystal,*) xf(jatoms),yf(jatoms),zf(jatoms)
        end do
      end do
    write(icrystal,'(i0)') 0
    write(icrystal,'(a)') trim(config_title)

    close(icrystal)
    return

  end subroutine write_crystal

!===============================================================================


  subroutine write_www
! ====================

    use structure_data
    use arguments
    use annotations
    use channel_numbers

    implicit none

    integer :: i,j,n,nl,jatom
    character(len=132) :: buffer,bline
    double precision :: dx,dy,dz,dxo,dyo,dzo,r,box,rmax

    n = index(fileroot,'.')
    if (n>0) then
      open(iwww,file=trim(wfolder)//fileroot(1:n)//'www',form='formatted',status='unknown')
    else
      open(iwww,file=trim(wfolder)//trim(fileroot)//'.www',form='formatted',status='unknown')
    end if

    box = (a + b + c)/6.0d0
    write(buffer,'(f12.5)') box
    buffer = adjustl(buffer)
    buffer = 'BOX='//trim(buffer)
    write(iwww,'(a)') trim(buffer)
    write(iwww,'(a)') 'energy=0.0'
    do i = 1,natoms
      write(iwww,'(i0,3f9.3)') i-1,2.0d0*(xf(i)-0.5d0)*box,2.0d0*(yf(i)-0.5d0)*box, &
                                   2.0d0*(zf(i)-0.5d0)*box  
    end do
    
    rmax = 2.0d0*box*(sqrt(3.0d0)/4.0d0)/(real(natoms)/8.0d0)**(1.0d0/3.0d0)
    rmax = 1.1d0*rmax

    do i = 1,natoms
      n = neighbours(i,0)
      write(bline,'(i0)') i - 1
      do j = 1,n
         jatom = neighbours(i,j)
         dx = xf(i) - xf(jatom) + 1.5d0
         dx = dx - aint(dx) - 0.5d0
         dy = yf(i) - yf(jatom) + 1.5d0
         dy = dy - aint(dy) - 0.5d0
         dz = zf(i) - zf(jatom) + 1.5d0
         dz = dz - aint(dz) - 0.5d0
!         dxo = cell(1,1)*dx + cell(1,2)*dy + cell(1,3)*dz
!         dyo = cell(2,1)*dx + cell(2,2)*dy + cell(2,3)*dz
!         dzo = cell(3,1)*dx + cell(3,2)*dy + cell(3,3)*dz
         r = 2.0d0*box*sqrt(dx**2 + dy**2 + dz**2)
         if (r>rmax) cycle
         nl = len_trim(bline) + 3
         write(bline(nl:),'(i0)') jatom-1
      end do
      write(iwww,'(a)') trim(bline)
    end do

    close(iwww)
    return

  end subroutine write_www


!===============================================================================


  subroutine write_xtl
! ====================

!--------------------------------------------------------------------------
!
!     This subroutine writes the configuration in the XTL format
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use annotations
    use channel_numbers

    implicit none

    integer :: i,n

    n = index(fileroot,'.')
    if (n>0) then
      open(ixtl,file=trim(wfolder)//fileroot(1:n)//'xtl',form='formatted',status='unknown')
    else
      open(ixtl,file=trim(wfolder)//trim(fileroot)//'.xtl',form='formatted',status='unknown')
    endif

    write(ixtl,'(a)') 'TITLE '//trim(metadata_title)
    write(ixtl,'(a)') 'CELL'
    write(ixtl,'(6f12.6)') a,b,c,alpha,beta,gamma
    write(ixtl,'(a)') 'ATOMS'
    write(ixtl,'(a)') 'NAME X Y Z'
    do i = 1,natoms
      write(ixtl,'(a2,2x,3f10.5)') element(atom_type(i)),xf(i),yf(i),zf(i)
    end do
    write(ixtl,'(a)') 'EOF'
    close(ixtl)

    return

  end subroutine write_xtl


!===============================================================================

  subroutine write_crush
! ======================

    use structure_data
    use arguments
    use annotations
    use channel_numbers

    implicit none

    integer :: i,j,n,iel,nrus,icount,ilink
    integer, allocatable :: ielements(:)

    double precision, allocatable :: bondmax(:)

    character(len=132) :: buffer,temp
    character(len=132), allocatable :: ru_atoms(:)
    character(len=2) :: link_element

! Sort numbers of atoms of each type
    ntypes = 0
    do iel = -1,nelements
    if (n_elements(iel)>0) then
      ntypes = ntypes + 1
    end if
    end do
    if (.not.allocated(numoftype)) allocate(numoftype(ntypes))
    ntypes = 0
    do iel = -1,nelements
    if (n_elements(iel)>0) then
      ntypes = ntypes + 1
      numoftype(ntypes) = n_elements(iel)
    write(6,'(i0,2a)') n_elements(iel),' atoms of element ',trim(element(iel))
    end if
    end do

    n = index(fileroot,'.')
    if (n>0) then
      open(icrush,file=trim(wfolder)//fileroot(1:n)//'crsh',form='formatted',status='unknown')
    else
      open(icrush,file=trim(wfolder)//trim(fileroot)//'.crsh',form='formatted',status='unknown')
    end if

! Write title
    if (.not.lm_title) then
     write(6,'(a)', advance="no") 'Please give title for configuration: '
     read(5,'(a)') config_title
    end if
    write(icrush,'(a)') trim(adjustl(config_title))

! Write lattice type
    write(icrush,'(a,a)') centre,'    ! Lattice centre'

    write(icrush,'(i0,a)') ntypes,'    ! Number of different atom types'

    do iel = -1,nelements
    if (n_elements(iel)>0) then
      write(icrush,'(f9.4,a,i0)') element_mass(iel), &
         '    ! Mass of atom type ',iel
    end if
    end do

    write(icrush,'(a)') 'F   ! Do not average the inertia tensor of rigid units'
    write(6,'(a)', advance="no") 'Produce file for ANALYSE code (T/F, default = F): '
    read(5,*) buffer
    temp = adjustl(buffer)
    if ((temp(1:1)/='T').and.(temp(1:1)/='t')) temp(1:2) = 'F '
    if (temp(1:1)=='t') temp(1:1) ='T'
    write(6,'(a)', advance="no") 'Produce file for LOCAL code (T/F, default = F): '
    read(5,*) buffer
    temp(3:) = adjustl(buffer)
    if ((temp(3:3)/='T').and.(temp(3:3)/='t')) temp(3:3) = 'F'
    if (temp(3:3)=='t') temp(3:3) ='T'
    temp(4:4) = ' '
    if (temp(3:3)=='T') then
    write(6,'(a)', advance="no") 'Number of modes to be analysed using LOCALRUM: '
    read(5,*) n
    write(temp(5:),'(i0)') n
    end if
    write(icrush,'(a)') trim(temp)//'    ! Files for ANALYSE and LOCAL RUM, # eigenvectors'
    write(6,'(a)', advance="no") 'Number of eigenvectors printed in output file: '
    read(5,*) n
    write(icrush,'(i0,a)') n,'    ! Number of eigenvectors printed in output file'
    write(icrush,'(a)') '0 F   ! Print out all frequencies, in THz'
    write(icrush,'(a)') 'F F   ! For eigenvectors, whether to use phase angles and remove origin'
    write(icrush,'(a)') 'T     ! Calculate participation ratio'
    write(icrush,'(a)') '0     ! Number of neutron Brillouin zones'
    write(icrush,'(a)') '0 0 0 0 0 0'
    write(*,'(a)', advance="no") 'Number of types of rigid unit you want to consider: '
    read(5,*) nrus
    allocate(ru_atoms(nrus))
    allocate(bondmax(nrus))
    do i = 1,nrus
    write(6,'(a,i0,a)', advance="no") 'Elements in rigid unit ',i,': '
    read(5,'(a)') ru_atoms(i)
    write(6,'(a,i0,a)', advance="no") 'Maximum bond length in rigid unit ',i,': '
    read(5,*) bondmax(i)
    end do
    write(6,'(a)', advance="no") 'Link element: '
    read(5,*) link_element
    link_element = adjustl(link_element)
    ilink = -10
    do iel = -1,111
    if ((link_element==element(iel)).or.(link_element==elementuc(iel)).or.(link_element==elementlc(iel))) then
      ilink = iel
      exit
    end if
    end do
    if (ilink==-10) then
    write(6,'(a,a,a)') 'Element ',trim(buffer(1:2)),' not recognised ...'
    write(6,'(a)', advance="no") 'Please give its atomic number: '
    read(5,*) ilink
    if (ilink>110) stop 'Atom is not recognised'
    if (ilink<-1) stop 'Atom is not recognised'
    end if

! Form the rigid units here ...
   do i = 1,nrus
    buffer = trim(adjustl(ru_atoms(i)))
    icount = 0    ! Count the number of element types
    do while (len(trim(buffer))>0)
      icount = icount + 1
      buffer = trim(adjustl(buffer(3:)))
    end do
    write(6,'(a,i0,a,i0)') 'You have ',icount,' elements in unit ',i
    buffer = trim(adjustl(ru_atoms(i)))
    allocate(ielements(icount))
    do j = 1,icount
      ielements(j) = -10
      do iel = -1,111
        if ((buffer(1:2)==element(iel)).or.(buffer(1:2)==elementuc(iel)).or.(buffer(1:2)==elementlc(iel))) then
          ielements(j) = iel
          exit
        end if
      end do
      if (ielements(j)==-10) then
        write(6,'(a,a,a)') 'Element ',trim(buffer(1:2)),' not recognised ...'
        write(6,'(a)', advance="no") 'Please give its atomic number: '
        read(5,*) ielements(j)
        if (ielements(j)>110) stop 'Atom is not recognised'
        if (ielements(j)<-1) stop 'Atom is not recognised'
      end if
      buffer = trim(adjustl(buffer(3:)))
    end do
    if (allocated(ielements)) deallocate(ielements)
    end do
    close(icrush)
    return

  end subroutine write_crush

!===============================================================================

  subroutine write_distance_windows_file
! ======================================

!--------------------------------------------------------------------------
!
!     This subroutine writes the Distance Windows neighours file
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use annotations
    use utilities
    use channel_numbers

    implicit none

    integer :: i,j,iconstraint,ierror,iline,nlines,nconstraints,il,jatom,n, &
               nb,nc,nlinks
    
    !molnum,i,j,k,l,m,ifield,ifieldout,ierror,ibond,iangle,irigid,nbonds, &
    !           nangles,iline,nlines,n,nc,il,nlinks,nb,na,ntorsions,itorsion,ninversions, &
    !           iinversion,nrigids,iatom,jatom,katom,nrigid_bodies, &
    !           junk,ij,irbatom,latom,iconstraint,nconstraints,item
    integer :: nlabel(4),ipad(1)
!    integer, allocatable :: nrigid(:),ntemp(:),rigidtemp(:),nrigid_body(:), &
!                            natoms_in_rigid_body(:),jatoms(:)
    logical :: lcontinue,lokay
!    logical, allocatable :: lrigid(:)
!    character(len=1000),allocatable :: crigid(:),ctemp(:)
!    character(len=1000) :: cpad(1)
    character(len=800)  :: text,buffer,constraintinfo,constraintstuff
!    character(len=4000) :: angleinfo,anglestuff,torsioninfo,torsionstuff,inversioninfo, &
!                           inversionstuff
    character(len=4000),allocatable :: constraint_data(:),constraint_info(:,:)
    character(len=2) :: catom,celement(4)
    character(len=10) :: clabel(4)
    character(len=4) :: constraint_pair
    double precision :: rmin,rmax,dx,dy,dz,dxo,dyo,dzo,r
    double precision,allocatable :: rbmin(:),rbmax(:)

    open(idw,file=trim(dwfile),form='formatted',status='old')

! Count number of lines in fieldfile
    ierror = 0 ; nlines = 0
    do while (ierror==0)
    read(idw,'(a)',iostat=ierror) text
    if (ierror/=0) exit
    nlines = nlines + 1
    end do
    rewind(idw)

! Count number of constraints
    ierror = 0 ; nconstraints = 0
    do iline = 1,nlines
    read(idw,'(a)',iostat=ierror) text
    text = adjustl(text)
    if (index(text,'dw')==0) cycle
    nconstraints = nconstraints + 1
    text = adjustl(text(3:))
    n = index(text,' ') - 1   ! Locate the first element
    call extract_atom_name_label(text(1:n),celement(1),clabel(1))
    text = adjustl(text(n+1:))
    n = index(text,' ') - 1   ! Locate the second element
    call extract_atom_name_label(text(1:n),celement(2),clabel(2))
    if (celement(1)/=celement(2)) nconstraints = nconstraints + 1
    end do
    rewind(idw)

! Now obtain constraint information
    ierror = 0 ; iconstraint = 0
    if (nconstraints>0) then
    allocate(constraint_data(nconstraints),constraint_info(natoms,nconstraints))
    constraint_info = ''
    do iline = 1,nlines
      read(idw,'(a)',iostat=ierror) text
      if (index(text,'dw')==0) cycle
      buffer = adjustl(text(3:))
      iconstraint = iconstraint + 1
      constraint_data(iconstraint) = trim(buffer)
      n = index(buffer,' ') - 1   ! Locate the first element
      call extract_atom_name_label(buffer(1:n),celement(1),clabel(1))
      buffer = adjustl(buffer(n+1:))
      n = index(buffer,' ') - 1   ! Locate the first element
      call extract_atom_name_label(buffer(1:n),celement(2),clabel(2))
      if (celement(1)==celement(2)) cycle
      buffer = adjustl(buffer(n+1:))
      iconstraint = iconstraint + 1
      buffer = trim(celement(2))//' '//trim(celement(1))//' '//trim(adjustl(buffer))
      constraint_data(iconstraint) = trim(buffer)
    end do
    rewind(idw)
    end if

!	do i = 1,nconstraints
!	  write(6,'(a)') trim(constraint_data(i))
!	end do
!	stop

! Locate constraints
  do iconstraint = 1,nconstraints
    buffer = adjustl(constraint_data(iconstraint))
!    buffer = adjustl(buffer(3:))   ! After the constraint word
    n = index(buffer,' ') - 1   ! Locate the first element
    call extract_atom_name_label(buffer(1:n),celement(1),clabel(1))
    if (trim(clabel(1))=='*') then
      nlabel(1) = -2
    else if (len_trim(clabel(1))>0) then
      read(clabel(1),*,iostat=ierror) nlabel(1)
    else
      nlabel(1) = 0
    end if
    buffer = adjustl(buffer(n+1:))
    n = index(buffer,' ') - 1   ! Locate the first element
    call extract_atom_name_label(buffer(1:n),celement(2),clabel(2))
    if (trim(clabel(2))=='*') then
      nlabel(2) = -2
    else if (len_trim(clabel(2))>0) then
      read(clabel(2),*,iostat=ierror) nlabel(2)
    else
      nlabel(2) = 0
    end if
    buffer = adjustl(buffer(n+1:))
    read(buffer,*) rmin,rmax
    do i = 1,natoms
      if (trim(element(atom_type(i)))/=trim(celement(1))) cycle
      if ((atom_label(i)/=nlabel(1)).and.(trim(clabel(1))/='*')) cycle
      do jatom = 1,neighbours(i,0)
        j = neighbours(i,jatom)
        if (i==j) cycle
        if (element(atom_type(j))/=celement(2)) cycle
        if ((atom_label(j)/=nlabel(2)).and.(trim(clabel(2))/='*')) cycle
        dx = xf(i) - xf(j) + 1.5d0
        dx = dx - aint(dx) - 0.5d0
        dy = yf(i) - yf(j) + 1.5d0
        dy = dy - aint(dy) - 0.5d0
        dz = zf(i) - zf(j) + 1.5d0
        dz = dz - aint(dz) - 0.5d0
        dxo = cell(1,1)*dx + cell(1,2)*dy + cell(1,3)*dz
        dyo = cell(2,1)*dx + cell(2,2)*dy + cell(2,3)*dz
        dzo = cell(3,1)*dx + cell(3,2)*dy + cell(3,3)*dz
        r = dsqrt(dxo**2 + dyo**2 + dzo**2)
        if (r>rmax) cycle
        if (r<rmin) cycle
        n = len_trim(constraint_info(i,iconstraint)) + 2
        write(constraint_info(i,iconstraint)(n:),'(i0)') j
        n = len_trim(constraint_info(i,iconstraint)) + 2
        constraint_info(i,iconstraint)(n:n+1) = element(atom_type(j))
        n = len_trim(constraint_info(i,iconstraint)) + 2
        constraint_info(i,iconstraint)(n:n+1) = ';'
      end do
    end do
  end do

! Count the number of constraints per atom
  do i = 1,natoms
    do iconstraint = 1,nconstraints
      n = len_trim(constraint_info(i,iconstraint))
      nc = 0
      do j = 1,n
        if(constraint_info(i,iconstraint)(j:j)==';') nc = nc+1
      end do
      write(constraint_info(i,iconstraint)(n+2:),'(i0)') nc
    end do
  end do

! Count the total number of constraints
  nb = 0
  do iconstraint = 1,nconstraints
    do i = 1,natoms
      constraintinfo = adjustl(constraint_info(i,iconstraint))
      n = index(trim(constraintinfo),';',.true.)
      if (n==0) cycle
      n = n + 1
      read(constraintinfo(n:),*) nlinks
      nb = nb + nlinks
!write(101,'(a)') element(atom_type(i))//': '//constraintinfo
    end do
  end do

    n = index(fileroot,'.')
    if (n>0) then
      open(idwout,file=trim(wfolder)//fileroot(1:n)//'neigh',form='formatted',status='unknown')
    else
      open(idwout,file=trim(wfolder)//trim(fileroot)//'.neigh',form='formatted',status='unknown')
    end if

! Write out the neighbour file
  write(idwout,'(i0)') nb
  do i = 1,natoms
  buffer = ''
  nc = 0
    do iconstraint = 1,nconstraints
      constraintinfo = adjustl(constraint_info(i,iconstraint))
      n = index(trim(constraintinfo),';',.true.)
      if (n==0) cycle
      n = n + 1
      read(constraintinfo(n:),*) nlinks
      n = index(constraint_data(iconstraint),':') + 1
      constraintstuff = adjustl(constraint_data(iconstraint)(n:))
      do il = 1,nlinks
        read(constraintinfo,*) j
        n = index(constraintinfo,';') + 1
        constraintinfo = adjustl(constraintinfo(n:))
!        write(idwout,'(i0,2x,i0,2x,a)') i,j,trim(constraintstuff)
        write(buffer(len_trim(buffer)+2:),'(i0)') j
        nc = nc + 1
      end do
    end do
    write(idwout,'(i0)') nc
    write(idwout,'(a)') trim(adjustl(buffer))
    buffer = ''
  end do
  
  return

  end subroutine write_distance_windows_file


!===============================================================================

  subroutine analysis
! ===================

!--------------------------------------------------------------------------
!
!     This subroutine gives some analysis of a configuration using methods
!     from other subroutines
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use annotations
    use channel_numbers

    implicit none

    integer :: i,j,ierror,ibond,nbonds, &
             iline,nlines,n,nc,nlinks
    integer :: j1,j2,jl1,jl2,jatom,nlist
    integer,allocatable :: nbhist(:,:),jatoms(:),originlist(:)
    character(len=800)  :: text,buffer,bondinfo
    character(len=4000),allocatable :: bond_data(:),bond_info(:,:)
    character(len=4) :: bond_pair,fangle
    double precision :: pi,rad2deg,rmin,rmax,dx,dy,dz,dxo,dyo,dzo,r,r1,r2
    double precision :: dxo1,dyo1,dzo1,dxo2,dyo2,dzo2,cangle,angle,xjunk
   character(len=2), allocatable :: originatom(:)
   logical,allocatable :: luse(:)

   allocate(luse(natoms))
   luse = .true.

   if (llistfile) then
     luse = .false.
     open(ilist,file=trim(listfile),status='old',form='formatted')
       read(ilist,'(a)') buffer
       if (buffer(1:21)=='Number of particles =') then  ! We have reduced size AtomEye file
         read(buffer(22:),*) nlist
         do i = 1,10  ;  read(ilist,*)  ; end do
         allocate(originlist(nlist),originatom(nlist))
         do i = 1,nlist
           read(ilist,'(a)') buffer
           buffer = adjustl(buffer)
           n = index(buffer,' ')
           buffer = (adjustl(buffer(n:)))
           originatom(i) = buffer(1:2)
           buffer = adjustl(buffer(3:))
           read(buffer,*) xjunk,xjunk,xjunk,originlist(i)
           luse(originlist(i)) = .true.
         end do
       end if
     close(ilist)
   end if

    text = ''
    pi = dacos(-1.0d0) ; rad2deg = 180.0d0/pi
    if (lm_title) text = adjustl(trim(config_title))
    if (lm_material) text = trim(text)//'; material = '//adjustl(trim(metadata_material))
    if (len(trim(text))==0) text = 'Date = '//adjustl(trim(metadata_date))

! Open the analysis output file
    n = index(fileroot,'.')
    if (n>0) then
      open(ianal,file=trim(wfolder)//fileroot(1:n)//'anal',form='formatted',status='unknown')
    else
      open(ianal,file=trim(wfolder)//trim(fileroot)//'.anal',form='formatted',status='unknown')
    end if

    write(ianal,'(3f20.10)') cell(:,1)
    write(ianal,'(3f20.10)') cell(:,2)
    write(ianal,'(3f20.10)') cell(:,3)

! Open bond file
    open(ibondfile,file=trim(bondfile),form='formatted',status='old')
! Count number of lines in bondfile
    ierror = 0 ; nlines = 0
    do while (ierror==0)
    read(ibondfile,'(a)',iostat=ierror) text
    if (ierror/=0) exit
    nlines = nlines + 1
    end do
    rewind(ibondfile)

! Count number of bonds
    ierror = 0 ; nbonds = 0
    do iline = 1,nlines
    read(ibondfile,'(a)',iostat=ierror) text
    text = adjustl(text)
    if (index(text,'bond')>0) nbonds = nbonds + 1
    end do
    rewind(ibondfile)
! Now obtain bond information
    ierror = 0 ; ibond = 0
    if (nbonds>0) then
    allocate(bond_data(nbonds),bond_info(natoms,nbonds),nbhist(nbonds,0:20))
    nbhist = 0
    bond_info = ''
    do iline = 1,nlines
      read(ibondfile,'(a)',iostat=ierror) text
      text = adjustl(text)
      if (index(text,'bond')==0) cycle
      ibond = ibond + 1
      bond_data(ibond) = trim(text)
    end do
    end if

    close(ibondfile)

!    write(6,'(i0,a)') nrigids,' rigid term(s) read from bond information file'
    write(6,'(i0,a)') nbonds,' bond(s) read from bond information file'

    do ibond = 1,nbonds
    write(6,'(a)') trim(bond_data(ibond))
    end do

! Locate bonds
    do ibond = 1,nbonds
    buffer = adjustl(bond_data(ibond))
    buffer = adjustl(buffer(6:))   ! After the bond word
    bond_pair(1:2) = buffer(1:2)
    buffer = adjustl(buffer(3:))
    bond_pair(3:4) = buffer(1:2)
    buffer = adjustl(buffer(3:))
    read(buffer,*) rmin,rmax
    do i = 1,natoms
      if (.not.luse(i)) cycle
      if (lxyzlimits.and.(xf(i)<xyzlimits(1))) cycle
      if (lxyzlimits.and.(xf(i)>xyzlimits(2))) cycle
      if (lxyzlimits.and.(yf(i)<xyzlimits(3))) cycle
      if (lxyzlimits.and.(yf(i)>xyzlimits(4))) cycle
      if (lxyzlimits.and.(zf(i)<xyzlimits(5))) cycle
      if (lxyzlimits.and.(zf(i)>xyzlimits(6))) cycle
      if (element(atom_type(i))==bond_pair(1:2)) then
        do jatom = 1,neighbours(i,0)
          j = neighbours(i,jatom)
          if (i==j) cycle
          if (element(atom_type(j))==bond_pair(3:4)) then
!             if ( (element(atom_type(j))==element(atom_type(i))) .and. (i>j) ) cycle
             if (i==j) cycle
             dx = xf(i) - xf(j) + 1.5d0
             dx = dx - aint(dx) - 0.5d0
             dy = yf(i) - yf(j) + 1.5d0
             dy = dy - aint(dy) - 0.5d0
             dz = zf(i) - zf(j) + 1.5d0
             dz = dz - aint(dz) - 0.5d0
             dxo = cell(1,1)*dx + cell(1,2)*dy + cell(1,3)*dz
             dyo = cell(2,1)*dx + cell(2,2)*dy + cell(2,3)*dz
             dzo = cell(3,1)*dx + cell(3,2)*dy + cell(3,3)*dz
             r = dsqrt(dxo**2 + dyo**2 + dzo**2)
             if (r>rmax) cycle
             if (r<rmin) cycle
             n = len_trim(bond_info(i,ibond)) + 2
             write(bond_info(i,ibond)(n:),'(i0)') j
             n = len_trim(bond_info(i,ibond)) + 2
             bond_info(i,ibond)(n:n+1) = element(atom_type(j))
             n = len_trim(bond_info(i,ibond)) + 2
             bond_info(i,ibond)(n:n+1) = ';'
          end if
        end do
      end if
    end do
    end do

! Count the number of bonds per atom
    do i = 1,natoms
      if (.not.luse(i)) cycle
      if (lxyzlimits.and.(xf(i)<xyzlimits(1))) cycle
      if (lxyzlimits.and.(xf(i)>xyzlimits(2))) cycle
      if (lxyzlimits.and.(yf(i)<xyzlimits(3))) cycle
      if (lxyzlimits.and.(yf(i)>xyzlimits(4))) cycle
      if (lxyzlimits.and.(zf(i)<xyzlimits(5))) cycle
      if (lxyzlimits.and.(zf(i)>xyzlimits(6))) cycle
      do ibond = 1,nbonds
        n = len_trim(bond_info(i,ibond))
        nc = 0
        do j = 1,n
          if(bond_info(i,ibond)(j:j)==';') nc = nc+1
        end do
        write(bond_info(i,ibond)(n+2:),'(i0)') nc
      end do
    end do

! Analyse the number of bonds
    do ibond = 1,nbonds
    do i = 1,natoms
      if (.not.luse(i)) cycle
      if (lxyzlimits.and.(xf(i)<xyzlimits(1))) cycle
      if (lxyzlimits.and.(xf(i)>xyzlimits(2))) cycle
      if (lxyzlimits.and.(yf(i)<xyzlimits(3))) cycle
      if (lxyzlimits.and.(yf(i)>xyzlimits(4))) cycle
      if (lxyzlimits.and.(zf(i)<xyzlimits(5))) cycle
      if (lxyzlimits.and.(zf(i)>xyzlimits(6))) cycle
      bondinfo = adjustl(bond_info(i,ibond))
      n = index(trim(bondinfo),';',.true.)
      if (n==0) cycle
      n = n + 1
      read(bondinfo(n:),*) nlinks
      nbhist(ibond,nlinks) = nbhist(ibond,nlinks) + 1
    end do
    end do

! Print out numbers of bonds
   do i = 0,20
     write(ianal,*) i,(nbhist(ibond,i),ibond=1,nbonds)
   end do

! do i = 1,natoms
!   if (lxyzlimits.and.(xf(i)<xyzlimits(1))) cycle
!   if (lxyzlimits.and.(xf(i)>xyzlimits(2))) cycle
!   if (lxyzlimits.and.(yf(i)<xyzlimits(3))) cycle
!   if (lxyzlimits.and.(yf(i)>xyzlimits(4))) cycle
!   if (lxyzlimits.and.(zf(i)<xyzlimits(5))) cycle
!   if (lxyzlimits.and.(zf(i)>xyzlimits(6))) cycle
!   write(11,'(a)') trim(bond_info(i,1))
! end do

    close(ianal)

! Form the list of bond angles per bond type
   if (langles) then
     do ibond = 1,nbonds
       write(fangle,'(i0)') ibond
       n = index(fileroot,'.')
       if (n>0) then
         open(ianal,file=trim(wfolder)//fileroot(1:n)//'angles'//trim(fangle),form='formatted',status='unknown')
       else
         open(ianal,file=trim(wfolder)//trim(fileroot)//'.angles'//trim(fangle),form='formatted',status='unknown')
       end if
       do i = 1,natoms
         if (.not.luse(i)) cycle
         if (lxyzlimits.and.(xf(i)<xyzlimits(1))) cycle
         if (lxyzlimits.and.(xf(i)>xyzlimits(2))) cycle
         if (lxyzlimits.and.(yf(i)<xyzlimits(3))) cycle
         if (lxyzlimits.and.(yf(i)>xyzlimits(4))) cycle
         if (lxyzlimits.and.(zf(i)<xyzlimits(5))) cycle
         if (lxyzlimits.and.(zf(i)>xyzlimits(6))) cycle
         bondinfo = adjustl(bond_info(i,ibond))
         n = index(trim(bondinfo),';',.true.)
         if (n==0) cycle
         n = n + 1
         read(bondinfo(n:),*) nlinks
         allocate(jatoms(nlinks))
         do j = 1,nlinks
           read(bondinfo,*) jatoms(j)
           n = index(trim(bondinfo),';')
           bondinfo = bondinfo(n+1:)
         end do
         do jl1 = 1,nlinks-1
           j1 = jatoms(jl1)
           dx = xf(i) - xf(j1) + 1.5d0
           dx = dx - aint(dx) - 0.5d0
           dy = yf(i) - yf(j1) + 1.5d0
           dy = dy - aint(dy) - 0.5d0
           dz = zf(i) - zf(j1) + 1.5d0
           dz = dz - aint(dz) - 0.5d0
           dxo1 = cell(1,1)*dx + cell(1,2)*dy + cell(1,3)*dz
           dyo1 = cell(2,1)*dx + cell(2,2)*dy + cell(2,3)*dz
           dzo1 = cell(3,1)*dx + cell(3,2)*dy + cell(3,3)*dz
           r = dsqrt(dxo1**2 + dyo1**2 + dzo1**2)
           dxo1 = dxo1/r ; dyo1 = dyo1/r ; dzo1 = dzo1/r
           r1 = r
           do jl2 = jl1+1,nlinks
             j2 = jatoms(jl2)
             dx = xf(i) - xf(j2) + 1.5d0
             dx = dx - aint(dx) - 0.5d0
             dy = yf(i) - yf(j2) + 1.5d0
             dy = dy - aint(dy) - 0.5d0
             dz = zf(i) - zf(j2) + 1.5d0
             dz = dz - aint(dz) - 0.5d0
             dxo2 = cell(1,1)*dx + cell(1,2)*dy + cell(1,3)*dz
             dyo2 = cell(2,1)*dx + cell(2,2)*dy + cell(2,3)*dz
             dzo2 = cell(3,1)*dx + cell(3,2)*dy + cell(3,3)*dz
             r = dsqrt(dxo2**2 + dyo2**2 + dzo2**2)
             dxo2 = dxo2/r ; dyo2 = dyo2/r ; dzo2 = dzo2/r
             r2 = r
             cangle = dxo1*dxo2 + dyo1*dyo2 + dzo1*dzo2
             angle = acos(cangle)*rad2deg
             write(ianal,*) i,j1,j2,r1,r2,angle,cangle
           end do
         end do
         deallocate(jatoms)
       end do
     close(ianal)
     end do
   end if

    return

    if (allocated(bond_data)) deallocate(bond_data)
    if (allocated(bond_info)) deallocate(bond_info)
    if (allocated(nbhist)) deallocate(nbhist)

  end subroutine analysis


!===============================================================================


    subroutine dlpoly_analysis
!   ==========================

!--------------------------------------------------------------------------
!
!     This subroutine analyses a DLPOLY configuration file based on the
!     definitions given in the FIELD file
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use channel_numbers
    use utilities
    
    implicit none
    
    integer :: iatom,jatom,katom,matom
    integer :: i,n,ierror
    double precision :: pi,angle,cangle,dx,dy,dz,r,r1,r2,r3
    double precision :: u1,u2,ux1,uy1,uz1,ux2,uy2,uz2
    double precision :: dxo,dyo,dzo,dxo1,dyo1,dzo1,dxo2,dyo2,dzo2,dxo3,dyo3,dzo3
    character(len=132) :: buffer,buffer2
    logical :: lokay
    
    pi = acos(-1.0d0)
    
    inquire(file='FIELD',exist=lokay)
    if (.not.lokay) then
      write(6,'(a)') 'FIELD file not found, so no DL_POLY configuration analysis performed'
      return
    end if
    open(ifield,file='FIELD',form='formatted',status='old')
    
    ierror = 0
    field_read: do while (ierror==0)
      read(ifield,'(a)',iostat=ierror) buffer
      if (ierror/=0) exit field_read
      buffer = adjustl(buffer)
      if (len_trim(buffer)==0) cycle field_read
      n = index(buffer,' ') - 1
      select case (buffer(1:n))
      case ('RIGID')
        read(buffer(7:),*) n
        do i = 1,n
          read(ifield,*)
        end do
      case ('CONSTRAINTS')
        read(buffer(12:),*) n
        do i = 1,n
          read(ifield,*)
        end do
      case ('BONDS')
        open(ianal,file='DLanal.bonds',form='formatted',status='unknown')
        read(buffer(6:),*) n
        do i = 1,n
          read(ifield,'(a)') buffer2
          read(buffer2(5:),*) iatom,jatom
          dx = xf(iatom) - xf(jatom) + 1.5d0
          dx = dx - aint(dx) - 0.5d0
          dy = yf(iatom) - yf(jatom) + 1.5d0
          dy = dy - aint(dy) - 0.5d0
          dz = zf(iatom) - zf(jatom) + 1.5d0
          dz = dz - aint(dz) - 0.5d0
          dxo = cell(1,1)*dx + cell(1,2)*dy + cell(1,3)*dz
          dyo = cell(2,1)*dx + cell(2,2)*dy + cell(2,3)*dz
          dzo = cell(3,1)*dx + cell(3,2)*dy + cell(3,3)*dz
          r = dsqrt(dxo**2 + dyo**2 + dzo**2)
          write(ianal,'(2(i0,1x),f10.5)') iatom,jatom,r
        end do
        close(ianal)
      case ('ANGLES')
        read(buffer(7:),*) n
        open(ianal,file='DLanal.angles',form='formatted',status='unknown')
        do i = 1,n
          read(ifield,'(a)') buffer2
          read(buffer2(5:),*) jatom,iatom,katom
          dx = xf(iatom) - xf(jatom) + 1.5d0
          dx = dx - aint(dx) - 0.5d0
          dy = yf(iatom) - yf(jatom) + 1.5d0
          dy = dy - aint(dy) - 0.5d0
          dz = zf(iatom) - zf(jatom) + 1.5d0
          dz = dz - aint(dz) - 0.5d0
          dxo1 = cell(1,1)*dx + cell(1,2)*dy + cell(1,3)*dz
          dyo1 = cell(2,1)*dx + cell(2,2)*dy + cell(2,3)*dz
          dzo1 = cell(3,1)*dx + cell(3,2)*dy + cell(3,3)*dz
          r1 = dsqrt(dxo1**2 + dyo1**2 + dzo1**2)
          dx = xf(iatom) - xf(katom) + 1.5d0
          dx = dx - aint(dx) - 0.5d0
          dy = yf(iatom) - yf(katom) + 1.5d0
          dy = dy - aint(dy) - 0.5d0
          dz = zf(iatom) - zf(katom) + 1.5d0
          dz = dz - aint(dz) - 0.5d0
          dxo2 = cell(1,1)*dx + cell(1,2)*dy + cell(1,3)*dz
          dyo2 = cell(2,1)*dx + cell(2,2)*dy + cell(2,3)*dz
          dzo2 = cell(3,1)*dx + cell(3,2)*dy + cell(3,3)*dz
          r2 = dsqrt(dxo2**2 + dyo2**2 + dzo2**2)
          cangle = (dxo1*dxo2 + dyo1*dyo2 + dzo1*dzo2)/(r1*r2)
          angle = acos(cangle)*180.0d0/pi
          write(ianal,'(3(i0,1x),4f12.5)') iatom,jatom,katom,r1,r2,angle,cangle
        end do
        close(ianal)
      case ('DIHEDRALS')
        read(buffer(10:),*) n
        open(ianal,file='DLanal.dihedrals',form='formatted',status='unknown')
        do i = 1,n
          read(ifield,'(a)') buffer2
          read(buffer2(5:),*) iatom,jatom,katom,matom
          dx = xf(iatom) - xf(jatom) + 1.5d0
          dx = dx - aint(dx) - 0.5d0
          dy = yf(iatom) - yf(jatom) + 1.5d0
          dy = dy - aint(dy) - 0.5d0
          dz = zf(iatom) - zf(jatom) + 1.5d0
          dz = dz - aint(dz) - 0.5d0
          dxo1 = cell(1,1)*dx + cell(1,2)*dy + cell(1,3)*dz
          dyo1 = cell(2,1)*dx + cell(2,2)*dy + cell(2,3)*dz
          dzo1 = cell(3,1)*dx + cell(3,2)*dy + cell(3,3)*dz
          r1 = dsqrt(dxo1**2 + dyo1**2 + dzo1**2)
          dx = xf(jatom) - xf(katom) + 1.5d0
          dx = dx - aint(dx) - 0.5d0
          dy = yf(jatom) - yf(katom) + 1.5d0
          dy = dy - aint(dy) - 0.5d0
          dz = zf(jatom) - zf(katom) + 1.5d0
          dz = dz - aint(dz) - 0.5d0
          dxo2 = cell(1,1)*dx + cell(1,2)*dy + cell(1,3)*dz
          dyo2 = cell(2,1)*dx + cell(2,2)*dy + cell(2,3)*dz
          dzo2 = cell(3,1)*dx + cell(3,2)*dy + cell(3,3)*dz
          r2 = dsqrt(dxo2**2 + dyo2**2 + dzo2**2)
          dx = xf(katom) - xf(matom) + 1.5d0
          dx = dx - aint(dx) - 0.5d0
          dy = yf(katom) - yf(matom) + 1.5d0
          dy = dy - aint(dy) - 0.5d0
          dz = zf(katom) - zf(matom) + 1.5d0
          dz = dz - aint(dz) - 0.5d0
          dxo3 = cell(1,1)*dx + cell(1,2)*dy + cell(1,3)*dz
          dyo3 = cell(2,1)*dx + cell(2,2)*dy + cell(2,3)*dz
          dzo3 = cell(3,1)*dx + cell(3,2)*dy + cell(3,3)*dz
          r3 = dsqrt(dxo1**2 + dyo1**2 + dzo1**2)
          ux1 = dyo1*dzo2 - dyo2*dzo1
          uy1 = dzo1*dxo2 - dzo2*dxo1
          uz1 = dxo1*dyo2 - dxo2*dyo1
          u1 = sqrt(ux1**2 + uy1**2 + uz1**2)
          ux2 = dyo2*dzo3 - dyo3*dzo2
          uy2 = dzo2*dxo3 - dzo3*dxo2
          uz2 = dxo2*dyo3 - dxo3*dyo2
          u2 = sqrt(ux2**2 + uy2**2 + uz2**2)
          cangle = (ux1*ux2 + uy1*uy2 + uz1*uz2)/(u1*u2)
          angle = acos(cangle)*180.0d0/pi
          write(ianal,'(4(i0,1x),5f12.5)') iatom,jatom,katom,matom,r1,r2,r3,angle,cangle
        end do
        close(ianal)
      end select
    end do field_read
    
    close(ifield)
    close(ianal)
    
    return

    end subroutine dlpoly_analysis



!===============================================================================


    subroutine rmc_analysis
!   =======================

!--------------------------------------------------------------------------
!
!     This subroutine analyses an rmc6f configuration file based on the
!     definitions given in the bonds and triplets files
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use channel_numbers
    use utilities
    
    implicit none
    
    integer :: iatom,jatom,katom
    integer :: i,j,n,nbondlist,nlocal,nbondnumber,nbondtypes
    integer :: nangletypes,ntripletlist,ntripletnumber
    double precision :: pi,angle,cangle,dx,dy,dz,r,r1,r2,r3
    double precision :: u1,u2,ux1,uy1,uz1,ux2,uy2,uz2
    double precision :: dxo,dyo,dzo,dxo1,dyo1,dzo1,dxo2,dyo2,dzo2,dxo3,dyo3,dzo3
    character(len=132) :: buffer,buffer2
    character(len=132), allocatable :: bondsdata(:),tripletsdata(:)
    character(len=2) :: atom1, atom2, atom3
    logical :: lbonds,ltriplets
    logical, allocatable :: lbondused(:,:)
    
    pi = acos(-1.0d0)

    n = index(fileroot,'.') - 1
    inquire(file=trim(wfolder)//fileroot(1:n)//'.bonds',exist=lbonds)
    inquire(file=trim(wfolder)//fileroot(1:n)//'.triplets',exist=ltriplets)

    if ( (.not.lbonds) .and. (.not.ltriplets)) then
      write(6,'(a)') 'Neither a bonds nor a triples file exist for -rmcanal action'
      return
    end if
    if (.not.lbonds) then
      write(6,'(a)') 'No bonds file found, analysis only for the triplets file'
    end if
    if (.not.ltriplets) then
      write(6,'(a)') 'No triplets file found, analysis only for the bonds file'
    end if

!   Analysis of bonds file
    if (lbonds) then
      open(irbonds,file=trim(wfolder)//fileroot(1:n)//'.bonds',form='formatted',status='old')
      open(irbout,file=trim(wfolder)//fileroot(1:n)//'.rmcbonds',form='formatted',status='unknown')
      bonds_header: do
        read(irbonds,'(a)') buffer
        buffer = adjustl(buffer)
        if (buffer(1:10)=='..........') exit bonds_header
        n = index(buffer,'=')
        if (n==0) cycle bonds_header
        select case(trim(buffer(1:n-1)))
          case('Number of bonds')
            read(buffer(n+1:),*) nbondtypes
          case('Number of atoms')
            read(buffer(n+1:),*) nbondlist
        end select
      end do bonds_header
      allocate(bondsdata(nbondlist))
      do i = 1,nbondlist
        read(irbonds,'(a)') buffer
        buffer = adjustl(buffer)
        bondsdata(i) = trim(buffer)
      end do
      allocate(lbondused(nbondlist,nbondlist))
      lbonds = .false.
      bondlist: do i = 1,nbondlist
        buffer = bondsdata(i)
        n = index(buffer,';',back=.true.)
        read(buffer(n+1:),*) nlocal
        buffer = buffer(1:n-1)
        read(buffer,*) iatom
        n = index(buffer,' ')
        buffer = adjustl(buffer(n:))
        atom1 = buffer(1:2)
        buffer = adjustl(buffer(3:)) 
        read(buffer,*) nbondnumber
        n = index(buffer,'::')
        buffer2 = adjustl(buffer(n+2:))
        locallist: do j = 1,nlocal
          read(buffer2,*) jatom
          if (lbondused(iatom,jatom)) cycle locallist
          lbondused(iatom,jatom) = .true.
          lbondused(jatom,iatom) = .true.
          n = index(buffer2,' ')
          buffer2 = adjustl(buffer2(n:))
          atom2 = buffer2(1:2)
          n = index(buffer2,';')
          if (n>0) buffer2 = adjustl(buffer2(n+1:))
          dx = xf(iatom) - xf(jatom) + 1.5d0
          dx = dx - aint(dx) - 0.5d0
          dy = yf(iatom) - yf(jatom) + 1.5d0
          dy = dy - aint(dy) - 0.5d0
          dz = zf(iatom) - zf(jatom) + 1.5d0
          dz = dz - aint(dz) - 0.5d0
          dxo = cell(1,1)*dx + cell(1,2)*dy + cell(1,3)*dz
          dyo = cell(2,1)*dx + cell(2,2)*dy + cell(2,3)*dz
          dzo = cell(3,1)*dx + cell(3,2)*dy + cell(3,3)*dz
          r = dsqrt(dxo**2 + dyo**2 + dzo**2)
          write(irbout,'(2(a,1x),3(i0,1x),f10.5)') atom1,atom2,nbondnumber,iatom,jatom,r       
        end do locallist
      end do bondlist
      close(irbonds)
      close(irbout)
      deallocate(lbondused)
      deallocate(bondsdata)
    end if

    if (ltriplets) then
      open(irangles,file=trim(wfolder)//fileroot(1:n)//'.triplets',form='formatted',status='old')
      open(iraout,file=trim(wfolder)//fileroot(1:n)//'.rmcangles',form='formatted',status='unknown')
      triplets_header: do
        read(irangles,'(a)') buffer
        buffer = adjustl(buffer)
        if (buffer(1:10)=='..........') exit triplets_header
        n = index(buffer,'=')
        if (n==0) cycle triplets_header
        select case(trim(buffer(1:n-1)))
          case('Number of triplets')
            read(buffer(n+1:),*) nangletypes
          case('Number of atoms')
            read(buffer(n+1:),*) ntripletlist
        end select
      end do triplets_header
      allocate(tripletsdata(ntripletlist))
      do i = 1,ntripletlist
        read(irangles,'(a)') buffer
        read(irangles,'(a)')
        buffer = adjustl(buffer)
        tripletsdata(i) = trim(buffer)
      end do
      tripletlist: do i = 1,ntripletlist
        buffer = tripletsdata(i)
        n = index(buffer,';',back=.true.)
        if (n==0) cycle tripletlist
        read(buffer(n+1:),*) nlocal
        if (nlocal==0) cycle tripletlist
        buffer = buffer(1:n-1)
        read(buffer,*) iatom
        n = index(buffer,' ')
        buffer = adjustl(buffer(n:))
        atom1 = buffer(1:2)
        buffer = adjustl(buffer(3:)) 
        read(buffer,*) ntripletnumber
        n = index(buffer,'::')
        buffer2 = adjustl(buffer(n+2:))
        do j = 1,nlocal
          read(buffer2,*) jatom, atom2, katom, atom3
          n = index(buffer2,';')
          if (n>0) buffer2 = adjustl(buffer2(n+1:))
          dx = xf(iatom) - xf(jatom) + 1.5d0
          dx = dx - aint(dx) - 0.5d0
          dy = yf(iatom) - yf(jatom) + 1.5d0
          dy = dy - aint(dy) - 0.5d0
          dz = zf(iatom) - zf(jatom) + 1.5d0
          dz = dz - aint(dz) - 0.5d0
          dxo1 = cell(1,1)*dx + cell(1,2)*dy + cell(1,3)*dz
          dyo1 = cell(2,1)*dx + cell(2,2)*dy + cell(2,3)*dz
          dzo1 = cell(3,1)*dx + cell(3,2)*dy + cell(3,3)*dz
          r1 = dsqrt(dxo1**2 + dyo1**2 + dzo1**2)
          dx = xf(iatom) - xf(katom) + 1.5d0
          dx = dx - aint(dx) - 0.5d0
          dy = yf(iatom) - yf(katom) + 1.5d0
          dy = dy - aint(dy) - 0.5d0
          dz = zf(iatom) - zf(katom) + 1.5d0
          dz = dz - aint(dz) - 0.5d0
          dxo2 = cell(1,1)*dx + cell(1,2)*dy + cell(1,3)*dz
          dyo2 = cell(2,1)*dx + cell(2,2)*dy + cell(2,3)*dz
          dzo2 = cell(3,1)*dx + cell(3,2)*dy + cell(3,3)*dz
          r2 = dsqrt(dxo2**2 + dyo2**2 + dzo2**2)
          cangle = (dxo1*dxo2 + dyo1*dyo2 + dzo1*dzo2)/(r1*r2)
          angle = acos(cangle)*180.0d0/pi
          write(iraout,'(3(a,1x),4(i0,1x),4f12.5)') atom1,atom2,atom3,ntripletnumber,iatom,jatom,katom,r1,r2,angle,cangle
        end do
      end do tripletlist

      
    end if

    return

end subroutine rmc_analysis


!===============================================================================


    subroutine ylm
!   ==============

!--------------------------------------------------------------------------
!
!     This subroutine calculates the orthonormal orientational variables and
!     Kubic harmonics.
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use channel_numbers
    use utilities
    
    implicit none

    integer :: ierror, i, n, m, nrigid_bodies, j,iatom,imol, iatomstart
    character(200) :: text,buffer
    double precision, allocatable :: xcf(:), ycf(:), zcf(:)
    double precision :: dx,dy,dz,dxo,dyo,dzo,x,y,z,ktemp
    double precision :: rcell(3,3),volume,sk,qk,r,pi,value,error

    double precision :: y1(-1:1),y2(-2:2),y3(-3:3),y4(-4:4)
    double precision, allocatable :: y1sum(:,:),y2sum(:,:),y3sum(:,:),y4sum(:,:)
    double precision, allocatable :: kubic(:,:),kubicsq(:,:)
    integer, allocatable :: natoms_in_rigid_body(:),atom_in_rigid_body(:,:)
    character(len=4000),allocatable :: crigid(:)
    logical :: lexist
    
    pi = acos(-1.0d0)
    
!   Check if the molecule file exists and if not, call locate_molecules
    inquire(file='molecules.dat', exist = lexist)
    if (lexist) then
      continue
    else
      call locate_molecules
    end if 
    
    open(imolecule, file=trim(wfolder)//'molecules.dat')
    ierror = 0 
    
!   Read in the number of molecules   
    read(imolecule, '(a)', iostat = ierror) buffer
    if (index(buffer, 'Number of molecules =') /= 0) then
      n = index(buffer, '=') + 1
      read(buffer(n:), *) nrigid_bodies
    end if   
    
    if (allocated(natoms_in_rigid_body)) deallocate(natoms_in_rigid_body)
    allocate(natoms_in_rigid_body(nrigid_bodies))
    allocate(crigid(nrigid_bodies))
    allocate(xcf(nrigid_bodies), ycf(nrigid_bodies), zcf(nrigid_bodies))
    if (allocated(xco).and.allocated(yco).and.allocated(zco)) then
           continue
        else
           allocate (xco(nrigid_bodies),yco(nrigid_bodies),zco(nrigid_bodies))
    end if
   
!   Read in the atom number within each molecule and the coordination of centre of mass
    ierror = 0
    i = 0
    do while (ierror == 0) 
      read(imolecule,'(a)',iostat=ierror) buffer
      buffer = adjustl(buffer)
      if(ierror/=0) exit
      if(len_trim(buffer)==0) cycle
       
      if ((ichar(buffer(1:1))>=ichar('0')).and.(ichar(buffer(1:1))<=ichar('9'))) then
        i = i + 1
        read(buffer,*) natoms_in_rigid_body(i),xco(i),yco(i),zco(i)
        buffer = adjustl(buffer)
        n = index(buffer,' ') ; buffer = adjustl(buffer(n:)) !pass natoms_in_rigid_body(i) 
        n = index(buffer,' ') ; buffer = adjustl(buffer(n:)) !pass xco 
        n = index(buffer,' ') ; buffer = adjustl(buffer(n:)) !pass yco   
        n = index(buffer,' ') ; buffer = adjustl(buffer(n:)) !pass zco  
        crigid(i)=buffer !now crigid(i) contains the atom number for each atom in molecule(i)
      end if
    end do

    close(imolecule)
    
    allocate(atom_in_rigid_body(nrigid_bodies,maxval(natoms_in_rigid_body))) 
    
    do i = 1, nrigid_bodies
      read(crigid(i),*) atom_in_rigid_body(i,:) !now construct a matrix which contains all the atom number for all the molecules     
    end do   

     volume = cell(1,1)*cell(2,2)*cell(3,3) + cell(1,2)*cell(2,3)*cell(3,1) + &
				cell(1,3)*cell(2,1)*cell(3,2) - cell(1,1)*cell(2,3)*cell(3,2) - &
				cell(1,3)*cell(2,2)*cell(3,1) - cell(1,2)*cell(2,1)*cell(3,3)
  
    ! Calculate the Cartesian coordinates of centre of mass
	   rcell(1,1) = cell(2,2)*cell(3,3) - cell(2,3)*cell(3,2) ! correct
	   rcell(1,2) = cell(2,3)*cell(3,1) - cell(2,1)*cell(3,3) ! correct
	   rcell(1,3) = cell(2,1)*cell(3,2) - cell(2,2)*cell(3,1) ! correct
	   rcell(2,1) = cell(3,2)*cell(1,3) - cell(3,3)*cell(1,2) ! correct
	   rcell(2,2) = cell(3,3)*cell(1,1) - cell(3,1)*cell(1,3) ! correct
	   rcell(2,3) = cell(3,1)*cell(1,2) - cell(3,2)*cell(1,1) ! correct
	   rcell(3,1) = cell(1,2)*cell(2,3) - cell(1,3)*cell(2,2) ! correct
	   rcell(3,2) = cell(1,3)*cell(2,1) - cell(1,1)*cell(2,3) ! correct
	   rcell(3,3) = cell(1,1)*cell(2,2) - cell(1,2)*cell(2,1) ! correct

	   rcell = rcell/volume

     do i = 1, nrigid_bodies
       xcf(i) = xf(atom_in_rigid_body(i,1))
       ycf(i) = yf(atom_in_rigid_body(i,1))
       zcf(i) = zf(atom_in_rigid_body(i,1))         
     end do
     
     do i = 1,nrigid_bodies 
       read(crigid(i),*) (atom_in_rigid_body(i,j),j=1,natoms_in_rigid_body(i))
!       weight = element_mass(atom_type(atom_in_rigid_body(i,1)))
       do m= 2,natoms_in_rigid_body(i)
            j = atom_in_rigid_body(i,m)
            dx = xcf(i)/(m-1) - xf(j) + 1.5d0
            dx = dx - aint(dx) - 0.5d0
            dy = ycf(i)/(m-1) - yf(j) + 1.5d0
            dy = dy - aint(dy) - 0.5d0
            dz = zcf(i)/(m-1) - zf(j) + 1.5d0
            dz = dz - aint(dz) - 0.5d0
            xf(j) = xcf(i)/(m-1) - dx
            yf(j) = ycf(i)/(m-1) - dy
            zf(j) = zcf(i)/(m-1) - dz
            xcf(i) = xcf(i) + xf(j)
            ycf(i) = ycf(i) + yf(j)
            zcf(i) = zcf(i) + zf(j)
 !          xcf(i) = xcf(i)*weight + xf(j)*element_mass(atom_type(j))
 !          ycf(i) = ycf(i)*weight + yf(j)*element_mass(atom_type(j))
 !          zcf(i) = zcf(i)*weight + zf(j)*element_mass(atom_type(j))
 !          weight = weight + element_mass(atom_type(j))
 !          xcf(i) = xcf(i)/weight 
 !          ycf(i) = ycf(i)/weight 
 !          zcf(i) = zcf(i)/weight 
       end do
        xcf(i) = xcf(i)/dble(natoms_in_rigid_body(i))
        ycf(i) = ycf(i)/dble(natoms_in_rigid_body(i))
        zcf(i) = zcf(i)/dble(natoms_in_rigid_body(i))
     end do 
     
    allocate(y1sum(-1:1,nrigid_bodies))
    allocate(y2sum(-2:2,nrigid_bodies))
    allocate(y3sum(-3:3,nrigid_bodies))
    allocate(y4sum(-4:4,nrigid_bodies))
    allocate(kubic(14,nrigid_bodies))
    allocate(kubicsq(14,nrigid_bodies))
    y1sum = 0.0d0
    y2sum = 0.0d0
    y3sum = 0.0d0
    y4sum = 0.0d0
    kubic = 0.0d0

    if (lmolcentre) then
      iatomstart = 1
    else
      iatomstart = 2
    end if
    molecule_loop: do imol = 1,nrigid_bodies
      atom_loop: do iatom = iatomstart,natoms_in_rigid_body(imol)
        i = atom_in_rigid_body(imol,iatom)
        dx = xf(i) - xcf(imol) + 1.5d0
        dx = dx - aint(dx) - 0.5d0
        dy = yf(i) - ycf(imol) + 1.5d0
        dy = dy - aint(dy) - 0.5d0
        dz = zf(i) - zcf(imol) + 1.5d0
        dz = dz - aint(dz) - 0.5d0
        dxo = cell(1,1)*dx + cell(1,2)*dy + cell(1,3)*dz
        dyo = cell(2,1)*dx + cell(2,2)*dy + cell(2,3)*dz
        dzo = cell(3,1)*dx + cell(3,2)*dy + cell(3,3)*dz
        r = dsqrt(dxo**2 + dyo**2 + dzo**2)
        x = dxo/r
        y = dyo/r
        z = dzo/r

        y1(-1) = y*sqrt(0.75d0/pi)
        y1(0)  = z*sqrt(0.75d0/pi)
        y1(1)  = x*sqrt(0.75d0/pi)

        y2(-2) = x*y*sqrt(15.d0/pi)/2.0d0
        y2(-1) = y*z*sqrt(15.d0/pi)/2.0d0
        y2(0)  = (2.0d0*z**2 - x**2 - y**2)*sqrt(5.0d0/pi)/4.0d0
        y2(1)  = z*x*sqrt(15.d0/pi)/2.0d0
        y2(2)  = (x**2 - y**2)*sqrt(15.d0/pi)/4.0d0

        y3(-3) = y*(3.0d0*x**2-y**2)*sqrt(17.5d0/pi)/4.0d0
        y3(-2) = x*y*z*sqrt(105.0d0/pi)/2.0d0
        y3(-1) = y*(4.0d0*z**2-x**2-y**2)*sqrt(10.5d0/pi)/4.0d0
        y3(0)  = z*(2.0d0*z**2-3.0d0*(x**2+y**2))*sqrt(7.0d0/pi)/4.0d0
        y3(1)  = x*(4.0d0*z**2-x**2-y**2)*sqrt(10.5d0/pi)/4.0d0
        y3(2)  = z*(x**2-y**2)*sqrt(105.0d0/pi)/4.0d0
        y3(3)  = x*(x**2-3.0d0*y**2)*sqrt(17.5d0/pi)/4.0d0

        y4(-4) = x*y*(x**2-y**2)*0.75d0*sqrt(35.0d0/pi)
        y4(-3) = y*z*(3.0d0*x**2-y**2)*0.75d0*sqrt(17.5d0/pi)
        y4(-2) = x*y*(7.0d0*z**2-1.0d0)*0.75d0*sqrt(5.0d0/pi)
        y4(-1) = y*z*(7.0d0*z**2-3.0d0)*0.75d0*sqrt(2.5d0/pi)
        y4(0)  = (35.0d0*z**4-30.0d0*z**2+3.0d0)*3.0d0/(16.0d0*sqrt(pi))
        y4(1)  = x*z*(7.0d0*z**2-3.0d0)*0.75d0*sqrt(2.5d0/pi)
        y4(2)  = (x**2-y**2)*(7.0d0*z**2-1.0d0)*0.375d0*sqrt(5.0d0/pi)
        y4(3)  = x*z*(x**2-3.0d0*y**2)*0.75d0*sqrt(17.5d0/pi)
        y4(4)  = (x**4+y**4-6.0d0*(x*y)**2)*3.0d0*sqrt(35.0d0/pi)/16.0d0

        y1sum(:,imol) = y1sum(:,imol) + y1
        y2sum(:,imol) = y2sum(:,imol) + y2
        y3sum(:,imol) = y3sum(:,imol) + y3
        y4sum(:,imol) = y4sum(:,imol) + y4

! Note that what I call K(12) and K(13) are the pair K(12,2) and K(12,3)
        qk = x**4 + y**4 + z**4
        sk = (x*y*z)**2
        ktemp = sqrt(21.0d0) * (5.0d0*qk-3.0d0)/4.0d0
        kubic(4,imol) = kubic(4,imol) + ktemp
        kubicsq(4,imol) = kubic(4,imol) + ktemp**2
        ktemp = sqrt(6.5d0) * (462.d0*sk + 21.0d0*qk - 17.0d0)/8.0d0
        kubic(6,imol) = kubic(6,imol) + ktemp
        kubicsq(6,imol) = kubic(6,imol) + ktemp**2
        ktemp = sqrt(561.0d0) * (65.0d0*qk**2 - 208.0d0*sk - 94.0d0*qk + 33.0d0)/32.0d0
        kubic(8,imol) = kubic(8,imol) + ktemp
        kubicsq(8,imol) = kubic(8,imol) + ktemp**2
        ktemp = sqrt(227.5d0) * (7106.0d0*qk*sk + 187.0d0*qk**2 - 3190.0d0*sk &
                - 264.0d0*qk + 85.0d0)/64.0d0
        kubic(10,imol) = kubic(10,imol) + ktemp
        kubicsq(10,imol) = kubic(10,imol) + ktemp**2
        ktemp = sqrt(11.0d0/41.0d0)*(3.0d0/128.0d0)*(1025.0d0*qk**3 - 2704156.0d0*sk**2 &
            + 352716.0d0*qk*sk + 4199.0d0*qk**2 - 232492.0d0*sk - 6526.0d0*qk + 2423.0d0)
        kubic(12,imol) = kubic(12,imol) + ktemp
        kubicsq(12,imol) = kubic(12,imol) + ktemp**2
        ktemp = sqrt(676039.0d0/246.0d0)*(1.0d0/128.0d0) * (1025.0d0*qk**3 - 16212.0d0*sk**2 &
                  - 8532.0d0*qk*sk - 2298.0d0*qk**2 + 4884.0d0*sk + 1677.0d0*qk - 396.0d0)
        kubic(13,imol) = kubic(13,imol) + ktemp
        kubicsq(13,imol) = kubic(13,imol) + ktemp**2
        ktemp = sqrt(51765.0d0)*(15.0d0/256.0d0) * (1311.0d0*qk**2*sk + &
          (437.0d0/18.0d0)*qk**3 - (6992.0d0/3.0d0)*sk**2 - 1573.2d0*qk*sk - (1577.0d0/30.0d0)*qk**2 + &
          (1501.0d0/3.0d0)*sk + (1109.0d0/30.0d0)*qk - 8.5d0)
        kubic(14,imol) = kubic(14,imol) + ktemp
        kubicsq(14,imol) = kubic(14,imol) + ktemp**2

      end do atom_loop

      kubic(4,imol) = kubic(4,imol)/natoms_in_rigid_body(imol)
      kubicsq(4,imol) = kubicsq(4,imol)/natoms_in_rigid_body(imol)
      kubic(6,imol) = kubic(6,imol)/natoms_in_rigid_body(imol)
      kubicsq(6,imol) = kubicsq(6,imol)/natoms_in_rigid_body(imol)
      kubic(8,imol) = kubic(8,imol)/natoms_in_rigid_body(imol)
      kubicsq(8,imol) = kubicsq(8,imol)/natoms_in_rigid_body(imol)
      kubic(10,imol) = kubic(10,imol)/natoms_in_rigid_body(imol)
      kubicsq(10,imol) = kubicsq(10,imol)/natoms_in_rigid_body(imol)
      kubic(12,imol) = kubic(12,imol)/natoms_in_rigid_body(imol)
      kubicsq(12,imol) = kubicsq(12,imol)/natoms_in_rigid_body(imol)
      kubic(13,imol) = kubic(13,imol)/natoms_in_rigid_body(imol)
      kubicsq(13,imol) = kubicsq(13,imol)/natoms_in_rigid_body(imol)
      kubic(14,imol) = kubic(14,imol)/natoms_in_rigid_body(imol)
      kubicsq(14,imol) = kubicsq(14,imol)/natoms_in_rigid_body(imol)

    end do molecule_loop

    n = index(fileroot,'.')
    if (n>0) then
      open(iylm,file=trim(wfolder)//fileroot(1:n)//'y1',form='formatted',status='unknown')
    else
      open(iylm,file=trim(wfolder)//trim(fileroot)//'.y1',form='formatted',status='unknown')
    end if    
    do imol = 1,nrigid_bodies
      write(iylm,*) imol,y1sum(:,imol)
    end do
    close(iylm)

    n = index(fileroot,'.')
    if (n>0) then
      open(iylm,file=trim(wfolder)//fileroot(1:n)//'y2',form='formatted',status='unknown')
    else
      open(iylm,file=trim(wfolder)//trim(fileroot)//'.y2',form='formatted',status='unknown')
    end if    
    do imol = 1,nrigid_bodies
      write(iylm,*) imol,y2sum(:,imol)
    end do
    close(iylm)

    n = index(fileroot,'.')
    if (n>0) then
      open(iylm,file=trim(wfolder)//fileroot(1:n)//'y3',form='formatted',status='unknown')
    else
      open(iylm,file=trim(wfolder)//trim(fileroot)//'.y3',form='formatted',status='unknown')
    end if    
    do imol = 1,nrigid_bodies
      write(iylm,*) imol,y3sum(:,imol)
    end do
    close(iylm)

    n = index(fileroot,'.')
    if (n>0) then
      open(iylm,file=trim(wfolder)//fileroot(1:n)//'y4',form='formatted',status='unknown')
    else
      open(iylm,file=trim(wfolder)//trim(fileroot)//'.y4',form='formatted',status='unknown')
    end if    
    do imol = 1,nrigid_bodies
      write(iylm,*) imol,y4sum(:,imol)
    end do
    close(iylm)

    n = index(fileroot,'.')
    if (n>0) then
      open(iylm,file=trim(wfolder)//fileroot(1:n)//'ylm',form='formatted',status='unknown')
    else
      open(iylm,file=trim(wfolder)//trim(fileroot)//'.ylm',form='formatted',status='unknown')
    end if    
    do i = 4,14
      select case (i)
        case (4,6,8,10,14)
          value = sum(kubic(i,:))/nrigid_bodies
          error = sum(kubicsq(i,:))/nrigid_bodies
          error = sqrt((error-value**2)/nrigid_bodies)
          write(iylm,'(a,i0,a,f10.4,a,f10.4)') 'Kubic function ',i,': ',value,'+/-',error
        case (12)
          value = sum(kubic(i,:))/nrigid_bodies
          error = sum(kubicsq(i,:))/nrigid_bodies
          error = sqrt((error-value**2)/nrigid_bodies)
          write(iylm,'(a,f10.4,a,f10.4)') 'Kubic function 12,2: ',value,'+/-',error
        case (13)
          value = sum(kubic(i,:))/nrigid_bodies
          error = sum(kubicsq(i,:))/nrigid_bodies
          error = sqrt((error-value**2)/nrigid_bodies)
          write(iylm,'(a,f10.4,a,f10.4)') 'Kubic function 12,3: ',value,'+/-',error
      end select
    end do
    do i = -1,1
      write(iylm,'(a,i0,a,2f10.4)') 'Y1 and Y1^2 averages for m =  ',i,': ', &
            sum(y1sum(i,:))/nrigid_bodies,sum(y1sum(i,:)**2)/nrigid_bodies
    end do
    do i = -2,2
      write(iylm,'(a,i0,a,2f10.4)') 'Y2 and Y2^2 averages for m =  ',i,': ', &
            sum(y2sum(i,:))/nrigid_bodies,sum(y2sum(i,:)**2)/nrigid_bodies
    end do
    do i = -3,3
      write(iylm,'(a,i0,a,2f10.4)') 'Y3 and Y3^2 averages for m =  ',i,': ', &
            sum(y3sum(i,:))/nrigid_bodies,sum(y3sum(i,:)**2)/nrigid_bodies
    end do
    do i = -4,4
      write(iylm,'(a,i0,a,2f10.4)') 'Y4 and Y4^2 averages for m =  ',i,': ', &
            sum(y4sum(i,:))/nrigid_bodies,sum(y4sum(i,:)**2)/nrigid_bodies
    end do
    do i = 1,nrigid_bodies
      write(iylm,*) y1sum(:,i),y2sum(:,i),y3sum(:,i),y4sum(:,i),kubic(:,i)
    end do
    
    close(iylm)
    
 end subroutine ylm    

!===============================================================================


  subroutine compare_files
! ========================

!--------------------------------------------------------------------------
!
!     This subroutine reads in a comparison file and the outputs all
!     differences in positions
!
!     so far only for rmc6f files
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use annotations
    use channel_numbers

    implicit none

    double precision :: dx,dy,dz,dxo,dyo,dzo,xft,yft,zft,r
    integer :: ios,i,j,n,matoms
    logical :: lcell,lvectors
    character(len=132) :: buffer,buffert
    character(len=2) :: atom_namet

    ios = 0

    open(icmp,file=trim(comparefile),status='old',form='formatted')
    open(icmpout,file=trim(comparefile)//'.out',status='unknown',form='formatted')
    read_loop: do while (ios==0)
      read(icmp,'(a)',iostat=ios) buffer
      buffer = adjustl(buffer)
      n = index(buffer,':')
!!!      if (n==0) n = len(trim(buffer)) + 1
      if (n>0) then
        do i = 1,n-1                          ! this bit converts to uppercase
          if ((ichar(buffer(i:i))>=ichar('a')).and.(ichar(buffer(i:i))<=ichar('z'))) &
            buffer(i:i) = char(ichar(buffer(i:i)) - ichar('a') + ichar('A'))
        end do
      end if
      n = index(buffer,':')
      if (index(buffer,'NUMBER OF TYPES OF ATOMS')>0) then
        cycle read_loop
      end if
      if (index(buffer,'ATOM TYPES PRESENT')>0) then
        cycle read_loop
      end if
      if (index(buffer,'NUMBER OF EACH ATOM TYPE')>0) then
        cycle read_loop
      end if
      if (index(buffer,'NUMBER OF ATOMS')>0) then
        n = index(buffer,':') + 1
        read(buffer(n:),*) matoms
        if (matoms/=natoms) then
          write(6,'(a)') 'Comparison is impossible because number of atoms does not match'
          stop
        end if
        cycle read_loop
      end if
      if_atoms: if (index(buffer,'ATOMS')>0) then
! Note that we are reading lines that look like
!     1   Na  [1]  0.000000    0.000000    0.000000     1   0   0   0
! The first integer is optional, so we need to test for its existence.
! The integer in square brackets is also optional, and we test for the bracket.
! Moreover, the last four integers are also optional, and again we need to check they exist
! There are two blocks of text, one for the data file defining a value for natoms and one
! where we need to count the number of atoms. We use the technique of looking for leading white
! space and then removing it using the ADJUSTL function.
          atom_label = 0   ! assign for the case where no label is given
          atom_loop_1: do i = 1,natoms
            read(icmp,'(a)',iostat=ios) buffer
            if (ios/=0) then
              write(6,'(a)') 'Exit> End of rmc6f file encountered during atom read'
              stop
            end if
            buffer = adjustl(buffer)                  ! first check for atom number
            buffert = buffer
            if ( (ichar(buffer(1:1))>=ichar('0')).and.(ichar(buffer(1:1))<=ichar('9')) ) then
!!!              read(buffer,*) atom_number(i)           ! read atom number
              n = index(buffer,' ')                   ! move past the atom number
              buffer = adjustl(buffer(n:))
!!!            else
!!!              atom_number(i) = i                      ! otherwise assign atom number
            end if
            atom_namet = buffer(1:2)           ! read element label
            buffer = buffer(3:)                       ! move past the element symbol
!***
!***  Define element_symbol(:), atom_number(:), atom_label(:)
!***
            buffer = adjustl(buffer)
            if (buffer(1:1)=='[') then                ! this bit to read atom label
              n = index(buffer,']')
              buffer = buffer(n+1:)
            end if
            read(buffer,*) xft,yft,zft          ! read x,y,z
            
            dx = xf(i) - xft + 1.5d0
            dx = dx - aint(dx) - 0.5d0
            dy = yf(i) - yft + 1.5d0
            dy = dy - aint(dy) - 0.5d0
            dz = zf(i) - zft + 1.5d0
            dz = dz - aint(dz) - 0.5d0
            dxo = cell(1,1)*dx + cell(1,2)*dy + cell(1,3)*dz
            dyo = cell(2,1)*dx + cell(2,2)*dy + cell(2,3)*dz
            dzo = cell(3,1)*dx + cell(3,2)*dy + cell(3,3)*dz
            r = dsqrt(dxo**2 + dyo**2 + dzo**2)
            write(icmpout,'(i0,1x,a2,1x,a2,7f9.4)') i,element(atom_type(i)),atom_namet,r,xf(i),yf(i),zf(i),xft,yft,zft
          end do atom_loop_1
      end if if_atoms
     end do read_loop
     
     close(icmp)
     close(icmpout)

    return

  end subroutine compare_files


!===============================================================================


  subroutine align_files
! ======================

!--------------------------------------------------------------------------
!
!     This subroutine reads in a reference file and then looks for systematic
!     offsets that it tries to correct for. For example, if one fractional
!     coordinate is shifted by one unit cell along one axis as apparently
!     found by comparison with the reference file, then the coordinates will
!     be shifted.
!
!     so far only for rmc6f files
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use annotations
    use channel_numbers

    implicit none

    double precision :: dx,dy,dz,dxo,dyo,dzo,xft,yft,zft,r
    integer :: ios,i,j,n,matoms,ndx,ndy,ndz
    logical :: lcell,lvectors
    character(len=132) :: buffer,buffert
    character(len=2) :: atom_namet

    ios = 0

    open(icmp,file=trim(alignfile),status='old',form='formatted')
    open(icmpout,file=trim(alignfile)//'.out',status='unknown',form='formatted')
    read_loop: do while (ios==0)
      read(icmp,'(a)',iostat=ios) buffer
      buffer = adjustl(buffer)
      n = index(buffer,':')
!!!      if (n==0) n = len(trim(buffer)) + 1
      if (n>0) then
        do i = 1,n-1                          ! this bit converts to uppercase
          if ((ichar(buffer(i:i))>=ichar('a')).and.(ichar(buffer(i:i))<=ichar('z'))) &
            buffer(i:i) = char(ichar(buffer(i:i)) - ichar('a') + ichar('A'))
        end do
      end if
      n = index(buffer,':')
      if (index(buffer,'NUMBER OF TYPES OF ATOMS')>0) then
        cycle read_loop
      end if
      if (index(buffer,'ATOM TYPES PRESENT')>0) then
        cycle read_loop
      end if
      if (index(buffer,'NUMBER OF EACH ATOM TYPE')>0) then
        cycle read_loop
      end if
      if (index(buffer,'NUMBER OF ATOMS')>0) then
        n = index(buffer,':') + 1
        read(buffer(n:),*) matoms
        if (matoms/=natoms) then
          write(6,'(a)') 'Comparison is impossible because number of atoms does not match'
          stop
        end if
        cycle read_loop
      end if
      if_atoms: if (index(buffer,'ATOMS')>0) then
! Note that we are reading lines that look like
!     1   Na  [1]  0.000000    0.000000    0.000000     1   0   0   0
! The first integer is optional, so we need to test for its existence.
! The integer in square brackets is also optional, and we test for the bracket.
! Moreover, the last four integers are also optional, and again we need to check they exist
! There are two blocks of text, one for the data file defining a value for natoms and one
! where we need to count the number of atoms. We use the technique of looking for leading white
! space and then removing it using the ADJUSTL function.
          atom_label = 0   ! assign for the case where no label is given
          atom_loop_1: do i = 1,natoms
            read(icmp,'(a)',iostat=ios) buffer
            if (ios/=0) then
              write(6,'(a)') 'Exit> End of rmc6f file encountered during atom read'
              stop
            end if
            buffer = adjustl(buffer)                  ! first check for atom number
            buffert = buffer
            if ( (ichar(buffer(1:1))>=ichar('0')).and.(ichar(buffer(1:1))<=ichar('9')) ) then
!!!              read(buffer,*) atom_number(i)           ! read atom number
              n = index(buffer,' ')                   ! move past the atom number
              buffer = adjustl(buffer(n:))
!!!            else
!!!              atom_number(i) = i                      ! otherwise assign atom number
            end if
            atom_namet = buffer(1:2)           ! read element label
            buffer = buffer(3:)                       ! move past the element symbol
!***
!***  Define element_symbol(:), atom_number(:), atom_label(:)
!***
            buffer = adjustl(buffer)
            if (buffer(1:1)=='[') then                ! this bit to read atom label
              n = index(buffer,']')
              buffer = buffer(n+1:)
            end if
            read(buffer,*) xft,yft,zft          ! read x,y,z
            
            dx = xf(i) - xft + 1.5d0
            dx = dx - aint(dx) - 0.5d0
            dx = dx*ncell(1)
            ndx = nint(dx)
            dy = yf(i) - yft + 1.5d0
            dy = dy - aint(dy) - 0.5d0
            dy = dy*ncell(2)
            ndy = nint(dy)
            dz = zf(i) - zft + 1.5d0
            dz = dz - aint(dz) - 0.5d0
            dz = dz*ncell(3)
            ndz = nint(dz)
            if (ndx/=0) xf(i) = xf(i) - dble(ndx)/dble(ncell(1))
            if (ndy/=0) yf(i) = yf(i) - dble(ndy)/dble(ncell(2))
            if (ndz/=0) zf(i) = zf(i) - dble(ndz)/dble(ncell(3))
!write(601,'(4(i0,1x),6f8.4)') i,ndx,ndy,ndz,xf(i),xf(i),zf(i),xft,yft,zft
          end do atom_loop_1
      end if if_atoms
     end do read_loop
     
     close(icmp)

    return

  end subroutine align_files

  
!===============================================================================


  subroutine crystalmakerbonds
! ============================

! Thus subroutine takes a bondlist file generated by CrystalMaker from the same
! configuration, and prints a list of bond lengths. The intention is that the
! user can generate a list of bonds from an ideal configuration and then generate
! a set of distances corresponding the the same ordered list.

    use structure_data
    use arguments
    use annotations
    use utilities
    use channel_numbers

    implicit none

    integer :: i,j,n,n1,n2,nbonds,ierror
    integer,allocatable :: iatom(:),jatom(:),ndipole(:)
    double precision :: dx,dy,dz,dxo,dyo,dzo,r
    double precision :: thetax,thetay,thetaz,phix,phiy,phiz,rad2deg
    double precision, allocatable :: dipole(:,:),distance(:)
    character(len=80) :: buffer
    logical, allocatable :: lrelevant(:)
    
    allocate(dipole(natoms,3))
    allocate(lrelevant(natoms))
    allocate(ndipole(natoms))
    
    rad2deg = 180.0d0/acos(-1.0d0)

!   Open the CrystalMaker bonds file and first count the number of bonds
    open(icmb,file=trim(cmbondfile),form='formatted',status='old')
    do i = 1,6
      read(icmb,*)
    end do
    ierror = 0
    nbonds = 0
    do while (ierror==0)
      read(icmb,*,iostat=ierror) i
      if (ierror/=0) exit
      nbonds = nbonds + 1
    end do
    allocate(iatom(nbonds),jatom(nbonds))
    allocate(distance(nbonds))
    rewind(icmb)
!   Now extract the atom numbers
    do i = 1,6
      read(icmb,*)
    end do
    do i = 1,nbonds
      read(icmb,'(a)') buffer
      buffer = adjustl(buffer)
      n = index(buffer,' ')
      buffer = adjustl(buffer(n:))  ! remove first number
      n1 = index(buffer,' ') - 1   ! Find length of first element name
      n = index(buffer,'-') + 1
      buffer = adjustl(buffer(n:))  ! remove first element name and dash
      n2 = index(buffer,' ') - 1   ! Find length of second element name
      buffer = adjustl(buffer(n2+1:))  ! remove second element name
      read(buffer(n1+1:),*) iatom(i)   ! read index of first element
      n = index(buffer,' ')
      buffer = adjustl(buffer(n:))  ! remove first element label
      read(buffer(n2+1:),*) jatom(i)   ! read index of first element
      n = index(buffer,' ')
      buffer = adjustl(buffer(n:))  ! leave the bond distance
      read(buffer,*) distance(i)
    end do
    close(icmb)
  
 !  Compute bond lengths and print to a file
    n = index(fileroot,'.')
    if (n>0) then
      open(icmb,file=trim(wfolder)//fileroot(1:n)//'bonds',form='formatted',status='unknown')
      open(itmp1,file=trim(wfolder)//fileroot(1:n)//'bonvec',form='formatted',status='unknown')
      open(itmp2,file=trim(wfolder)//fileroot(1:n)//'bonang',form='formatted',status='unknown')
    else
      open(icmb,file=trim(wfolder)//trim(fileroot)//'.bonds',form='formatted',status='unknown')
      open(itmp1,file=trim(wfolder)//trim(fileroot)//'.bonvec',form='formatted',status='unknown')
      open(itmp2,file=trim(wfolder)//trim(fileroot)//'.bonang',form='formatted',status='unknown')
    end if
    dipole = 0.0d0
    ndipole = 0
    do n = 1,nbonds
      i = iatom(n)
      lrelevant(i) = .true.
      j = jatom(n)
      dx = xf(i) - xf(j) + 1.5d0
      dx = dx - aint(dx) - 0.5d0
      dy = yf(i) - yf(j) + 1.5d0
      dy = dy - aint(dy) - 0.5d0
      dz = zf(i) - zf(j) + 1.5d0
      dz = dz - aint(dz) - 0.5d0
      dxo = cell(1,1)*dx + cell(1,2)*dy + cell(1,3)*dz
      dyo = cell(2,1)*dx + cell(2,2)*dy + cell(2,3)*dz
      dzo = cell(3,1)*dx + cell(3,2)*dy + cell(3,3)*dz
      r = dsqrt(dxo**2 + dyo**2 + dzo**2)
!       if (abs(r-distance(n))>1.0) then
!         write(6,'(a)') 'Something wrong with the crystalmaker analysis'
!         write(6,'(a)') 'With file '//trim(fileroot)
!         write(6,'(a,i0)') 'for bond number ',n
!         write(6,'(a,2x,f7.3)') 'CrystalMaker distance = ',distance(n)
!         write(6,'(a,2x,f7.3)') 'data2config distance =  ',r
!         write(6,'(a,i0,a,i0)') 'Atoms are number ',i,' and ',j
!         write(6,'(a)') 'The two sets of fractional coordinates are:'
!         write(6,'(3f8.4)') xf(i),yf(i),zf(i)
!         write(6,'(3f8.4)') xf(j),yf(j),zf(j)
!         write(6,'(a)') 'The calculated separations of fractional coordinates are:'
!         write(6,'(3f8.4)') dx,dy,dz
!         write(6,'(a)') 'The calculated separations of orthogonal coordinates are:'
!         write(6,'(3f8.3)') dxo,dyo,dzo
!         close(icmb)
!         close(itmp1)
!         close(itmp2)
!       end if
      write(icmb,'(4(i0,2x),2(a,2x),f10.5)') n,i,j,reference_number(i),element(atom_type(i)), &
                                               element(atom_type(j)),r
      dipole(i,1) = dipole(i,1) + dxo
      dipole(i,2) = dipole(i,2) + dyo
      dipole(i,3) = dipole(i,3) + dzo
      ndipole(i) = ndipole(i) + 1
      write(itmp1,'(5(i0,2x),2(a,2x),7f10.5)') n,i,j,reference_number(i),reference_number(j), &
                element(atom_type(i)),element(atom_type(j)),dxo/r,dyo/r,dzo/r,dxo,dyo,dzo,r
      dxo = dxo/r ; dyo = dyo/r ; dzo = dzo/r
      thetax = rad2deg*acos(dxo) ; thetay = rad2deg*acos(dyo) ; thetaz = rad2deg*acos(dzo)
      r = dsqrt(dxo**2 + dyo**2)
      phiz = rad2deg*acos(dxo/r)
      if (dyo<0) phiz = 360.0d0-phiz
      r = dsqrt(dxo**2 + dzo**2)
      phiy = rad2deg*acos(dzo/r)
      if (dxo<0) phiy = 360.0d0-phiy
      r = dsqrt(dyo**2 + dzo**2)
      phix = rad2deg*acos(dyo/r)
      if (dzo<0) phix = 360.0d0-phix
      write(itmp2,'(4(i0,2x),2(a,2x),6f10.5)') n,i,j,reference_number(i),element(atom_type(i)),element(atom_type(j)), &
            thetax,phix,thetay,phiy,thetaz,phiz
    end do 
    close(icmb)
    close(itmp1)
    close(itmp2)
    
    n = index(fileroot,'.')
    if (n>0) then
      open(icmb,file=trim(wfolder)//fileroot(1:n)//'dpl',form='formatted',status='unknown')
    else
      open(icmb,file=trim(wfolder)//trim(fileroot)//'.dpl',form='formatted',status='unknown')
    end if
    do i = 1,natoms
      if (lrelevant(i)) then
        write(icmb,'(3(i0,2x),3f12.5)') i,reference_number(i),ndipole(i),dipole(i,:)/dble(ndipole(i))
      end if
    end do

    deallocate(iatom,jatom,ndipole)
    deallocate(dipole)

    return

  end subroutine crystalmakerbonds
  

!===============================================================================

  subroutine locate_molecules
! ===========================

!--------------------------------------------------------------------------
!
!     This subroutine locates molecules in a configuration
!
!--------------------------------------------------------------------------

    use structure_data
    use arguments
    use annotations
    use utilities
    use channel_numbers

    implicit none

    integer :: molnum,i,j,k,l,m,ierror,irigid,nbonds,ir, &
               iline,nlines,n,nc,il,nlinks,nb,na, &
               nrigids,iatom,jatom,katom,nrigid_bodies, &
               junk,ij,irbatom,latom,item, &
               nrigid_bodies_temp, nrigid_bodies_total
    integer :: nlabel(4),ipad(1)
    integer, allocatable :: nrigid(:),ntemp(:),rigidtemp(:),nrigid_body(:), &
                            natoms_in_rigid_body(:),jatoms(:)
    logical :: lcontinue,lokay,ltemplate
    logical, allocatable :: lrigid(:)
    character(len=1000),allocatable :: crigid(:),ctemp(:)
    character(len=1000) :: cpad(1)
    character(len=800)  :: text,buffer
    character(len=4000),allocatable :: rigid_data(:)
    character(len=4),allocatable :: rigid_pair(:)
    character(len=2) :: catom,celement(4)
    character(len=10) :: clabel(4)
    character(len=4) :: bond_pair,cb1,cb2
    double precision :: pi,rbond,dr,rmin,rmax,dx,dy,dz,dxo,dyo,dzo, &
             r,r1,r2,r3,r4,r12,r13,r23,r24,r34,tmp,tmp1,tmp2, &
             r1min,r1max,r2min,r2max,r3min,r3max,r4min,r4max,r34min,r34max, &
             dx1,dx2,dx3,dx4,dx12,dx13,dx23,dx24,dx34, &
             dy1,dy2,dy3,dy4,dy12,dy13,dy23,dy24,dy34, &
             dz1,dz2,dz3,dz4,dz12,dz13,dz23,dz24,dz34, &
             dx1o,dx2o,dx3o,dx4o,dx12o,dx13o,dx23o,dx24o,dx34o, &
             dy1o,dy2o,dy3o,dy4o,dy12o,dy13o,dy23o,dy24o,dy34o, &
             dz1o,dz2o,dz3o,dz4o,dz12o,dz13o,dz23o,dz24o,dz34o, &
             angmean,ang1,ang2,dang,angmin,angmax,ang1min,ang1max,ang2min,ang2max, &
             angle,angle1,angle2,weight
    double precision,allocatable :: rbmin(:),rbmax(:)

    text = ''
    pi = dacos(-1.0d0)
    if (lm_title) text = adjustl(trim(config_title))
    if (lm_material) text = trim(text)//'; material = '//adjustl(trim(metadata_material))
    if (len(trim(text))==0) text = trim(text)//'; Date = '//adjustl(trim(metadata_date))

! Now we construct the FIELD file if required
    open(ifield,file=trim(moleculefile),form='formatted',status='old')

! Count number of lines in moleculefile
    ierror = 0 ; nlines = 0
    do while (ierror==0)
    read(ifield,'(a)',iostat=ierror) text
    if (ierror/=0) exit
    nlines = nlines + 1
    end do
    rewind(ifield)

! Count number of bond lines
    ierror = 0 ; nrigids = 0
    do iline = 1,nlines
      read(ifield,'(a)',iostat=ierror) text
      text = adjustl(text)
      if (index(text,'bond')>0) nrigids = nrigids + 1
    end do
    rewind(ifield)

! Now obtain rigid term information
    ierror = 0 ; irigid = 0
    if (nrigids>0) then
    allocate(rigid_data(nrigids))
    do iline = 1,nlines
      read(ifield,'(a)',iostat=ierror) text
      text = adjustl(text)
      if (index(text,'bond')==0) cycle
      irigid = irigid + 1
      rigid_data(irigid) = trim(text)
    end do
    rewind(ifield)
    end if

    close(ifield)

! Open the FIELDout file
   open(ifieldout,file=trim(wfolder)//'molecules.dat',form='formatted',status='unknown')


! To find molecules:
! 1. Loop over all atoms, exit loop if atom is already allocated to a rigid body
! 2. Start a new molecule, increment nrigid_bodies by one, and start a new text string
!    crigid in which the first number is the atom number. Also create an integer array
!    called natoms_in_rigid_body(nrigid_bodies) that contains the number 1 for this entry.
!    I will likely need to reallocate this array but I now know this is easy to do.
! 3. Now do a loop that ends when we add no more atoms to the rigid body. For this loop,
!    set n = natoms_in_rigid_body(nrigid_bodies) and loop over n. We can end the loop
!    when at the end n still equals natoms_in_rigid_body(nrigid_bodies).
! 4. Loop over the number of atoms within the text string. For each atom, loop over all its
!    neighbours, moving on if their logical lrigid is already true, and then loop over
!    the rigid bonds specifiers to see whether they are part of the rigid body. If they are
!    a) set lrigid to true; b) increment natoms_in_rigid_body(nrigid_bodies) by 1; c) add
!    the atom number to the end of the crigid text string.

     allocate(nrigid_body(natoms))
     allocate (lrigid(natoms))
     allocate (rigid_pair(nrigids),rbmin(nrigids),rbmax(nrigids))
     lrigid = .false.
     nrigid_body = 0
     nrigid_bodies = 0
     if (allocated(crigid)) crigid = ''
     cpad = ''
     do irigid = 1,nrigids
       buffer = adjustl(rigid_data(irigid))
       buffer = adjustl(buffer(6:))   ! After the rigid word
       bond_pair(1:2) = buffer(1:2)   ! Almost all this bit is the same as the bond part
       buffer = adjustl(buffer(3:))
       bond_pair(3:4) = buffer(1:2)
       rigid_pair(irigid) = bond_pair
       buffer = adjustl(buffer(3:))
       read(buffer,*) rbond,dr
       rbmin(irigid) = rbond - dr
       rbmax(irigid) = rbond + dr
     end do

     rigid_atom_loop_1: do iatom = 1,natoms
       if (lrigid(iatom)) cycle rigid_atom_loop_1
       nrigid_bodies = nrigid_bodies + 1   ! increment number of rigid bodies
       lrigid(iatom) = .true.
       if (allocated(nrigid)) then
         nrigid = reshape(nrigid,(/nrigid_bodies/),(/1/))
       else
         allocate(nrigid(1))
       end if
       if (allocated(natoms_in_rigid_body)) then
         natoms_in_rigid_body = reshape(natoms_in_rigid_body,(/nrigid_bodies/),(/1/))
         write(cpad(1),'(i0)') iatom
         crigid = reshape(crigid,(/nrigid_bodies/),cpad)
       else
         allocate(natoms_in_rigid_body(1))
         natoms_in_rigid_body(1) = 1
         allocate(crigid(1))
         write(crigid(1),'(i0)') iatom
       end if
       loop_over_list: do
         m = natoms_in_rigid_body(nrigid_bodies)
         if (allocated(jatoms)) deallocate(jatoms)
         allocate(jatoms(m))
         read(crigid(nrigid_bodies),*) jatoms
         loop_over_atoms_in_list: do item = 1,m
           i = jatoms(item)
           rigid_atom_loop_2: do jatom = 1,neighbours(i,0)
             j = neighbours(i,jatom)
             if (i==j) cycle rigid_atom_loop_2
             if (lrigid(j)) cycle rigid_atom_loop_2
             rigids_loop_2: do irigid = 1,nrigids
               cb1 = element(atom_type(i))//element(atom_type(j))
               cb2 = element(atom_type(j))//element(atom_type(i))
               if ((trim(cb1)/=trim(rigid_pair(irigid))).and.(trim(cb2)/=trim(rigid_pair(irigid)))) cycle rigids_loop_2
               dx = xf(i) - xf(j) + 1.5d0
               dx = dx - aint(dx) - 0.5d0
               dy = yf(i) - yf(j) + 1.5d0
               dy = dy - aint(dy) - 0.5d0
               dz = zf(i) - zf(j) + 1.5d0
               dz = dz - aint(dz) - 0.5d0
               dxo = cell(1,1)*dx + cell(1,2)*dy + cell(1,3)*dz
               dyo = cell(2,1)*dx + cell(2,2)*dy + cell(2,3)*dz
               dzo = cell(3,1)*dx + cell(3,2)*dy + cell(3,3)*dz
               r = dsqrt(dxo**2 + dyo**2 + dzo**2)
               if (r<rbmin(irigid)) cycle
               if (r>rbmax(irigid)) cycle
               natoms_in_rigid_body(nrigid_bodies) = natoms_in_rigid_body(nrigid_bodies) + 1
               n = len_trim(crigid(nrigid_bodies))+1
               write(crigid(nrigid_bodies)(n:),'(1x,i0)') j
               lrigid(j) = .true.
               lcontinue = .true.
             end do rigids_loop_2
           end do rigid_atom_loop_2
         end do loop_over_atoms_in_list
         if (m==natoms_in_rigid_body(nrigid_bodies)) exit loop_over_list
       end do loop_over_list
     end do rigid_atom_loop_1

! This bit exists because in the previous stuff we count individual atoms as
! rigid bodies, so here we need to decrement the total count
    nrigid_bodies_temp = nrigid_bodies
    do i = 1,nrigid_bodies
      if (natoms_in_rigid_body(i)==1) nrigid_bodies_temp = nrigid_bodies_temp - 1
    end do
    nrigid_bodies_total = nrigid_bodies
    nrigid_bodies = nrigid_bodies_temp

! Locate the centre of each molecule
     
     !guanqun:allocate the coordinates for molecule centres
      if (allocated(xco)) deallocate(xco)
      if (allocated(yco)) deallocate(yco)   
      if (allocated(zco)) deallocate(zco)
      allocate (xco(nrigid_bodies), yco(nrigid_bodies),zco(nrigid_bodies))
  ! Write out molecule data
    write(ifieldout,'(a,i0)') 'Number of molecules = ',nrigid_bodies
     !write molecule information including number of atoms in a molecule and their atom numbers
    ir = 0
     do i = 1,nrigid_bodies_total
      if (natoms_in_rigid_body(i)==1) cycle
      ir = ir + 1
      allocate(rigidtemp(natoms_in_rigid_body(i)))
      read(crigid(i),*) rigidtemp
      ! guanqun:calculate the molecule centre via averaging each atoms coordinates in a molecule
    !  xco(i)=0
     ! yco(i)=0
      !zco(i)=0
      !do m= 1, natoms_in_rigid_body(i)
      !xco(i)=xco(i)+xf(rigidtemp(m))
      !yco(i)=yco(i)+yf(rigidtemp(m))
      !zco(i)=zco(i)+zf(rigidtemp(m))
      !end do
      !m=natoms_in_rigid_body(i)
      !xco(i)=xco(i)/m
      !yco(i)=yco(i)/m
      !zco(i)=zco(i)/m
      
!      centre_of_mass = 0.0d0
!      weight = 0.0d0
      !guanqun:calculate the centre of mass
      xco(ir)=0.0d0
      yco(ir)=0.0d0
      zco(ir)=0.0d0
      weight= 0.0d0
      do m= 1, natoms_in_rigid_body(i)
      xco(ir)=xco(ir)+xo(rigidtemp(m))*element_mass(atom_type(rigidtemp(m)))
      yco(ir)=yco(ir)+yo(rigidtemp(m))*element_mass(atom_type(rigidtemp(m)))
      zco(ir)=zco(ir)+zo(rigidtemp(m))*element_mass(atom_type(rigidtemp(m)))
      weight = weight + element_mass(atom_type(rigidtemp(m)))
      end do
      xco(ir)=xco(ir)/weight
      yco(ir)=yco(ir)/weight
      zco(ir)=zco(ir)/weight
      
! Bootstrap this. Calculate the centre of first two atoms, then generate the nearest image of the third atom,
! find a new centre, and repeat for all atoms
! Then check that the centre is within the box still
! From the new centre generate the nearest atoms of the molecule and recalculate the xf and xo values
      write(ifieldout,'(i0,2x,3(f10.4,2x),16(i0,2x))') natoms_in_rigid_body(i), xco(ir), yco(ir), zco(ir), rigidtemp
      deallocate(rigidtemp)
    end do
  
  close(ifieldout)
  
  return
  
  end subroutine locate_molecules


!===============================================================================
  

  subroutine make_into_silica
! ===========================

!--------------------------------------------------------------------------
!
!     This subroutine expands a unit cell from a WWW configuration and adds
!     O in between
!
!--------------------------------------------------------------------------

    use arguments
    use annotations
    use structure_data
    use channel_numbers
    
    implicit none

    integer, allocatable :: atom_typet(:),atom_labelt(:),rnumt(:),rcellt(:,:)
    double precision,allocatable :: xft(:),yft(:),zft(:),xot(:),yot(:),zot(:)
    character(len=10), allocatable :: catom_labelt(:)
    character(len=4), allocatable :: atom_namet(:)
    
    integer :: i,j,k,jatom,neighbourlist(4),newatoms
    double precision :: dx,dy,dz,scale_factor
    character(len=80) :: buffer

    newatoms = natoms*3

    allocate(xft(natoms),yft(natoms),zft(natoms))
    allocate(xot(natoms),yot(natoms),zot(natoms))
    allocate(atom_namet(natoms))
    allocate(atom_typet(natoms))
    allocate(atom_labelt(natoms))
    allocate(rnumt(natoms),rcellt(natoms,3))
    xft = xf(1:natoms)
    yft = yf(1:natoms)
    zft = zf(1:natoms)
    xot = xo(1:natoms)
    yot = yo(1:natoms)
    zot = zo(1:natoms)
    atom_namet = atom_name(1:natoms)
    atom_typet = atom_type(1:natoms)
!    if (luselabels) then
!      catom_labelt = catom_label(1:natoms)
!    else
!      atom_labelt = atom_label(1:natoms)
!    end if
    atom_labelt = atom_label(1:natoms)
!    rnumt = reference_number(1:natoms)
!    rcellt = reference_cell(1:natoms,:)
    if (allocated(xf)) deallocate(xf,yf,zf)
    if (allocated(xo)) deallocate(xo,yo,zo)
    if (allocated(atom_name)) deallocate(atom_name)
    if (allocated(atom_type)) deallocate(atom_type)
    if (allocated(atom_label)) deallocate(atom_label)
    if (allocated(catom_label)) deallocate(catom_label)
    if (allocated(reference_number)) deallocate(reference_number)
    if (allocated(reference_cell)) deallocate(reference_cell)
    allocate(xf(newatoms),yf(newatoms),zf(newatoms))
    allocate(xo(newatoms),yo(newatoms),zo(newatoms))
    allocate(atom_name(newatoms))
    allocate(atom_type(newatoms))
!    if (luselabels) then
!      allocate(catom_label(newatoms))
!    else
!      allocate(atom_label(newatoms))
!    end if
    allocate(atom_label(newatoms))
    allocate(reference_number(newatoms))
    allocate(reference_cell(newatoms,3))
    xf(1:natoms) = xft
    yf(1:natoms) = yft
    zf(1:natoms) = zft
    xo(1:natoms) = xot
    yo(1:natoms) = yot
    zo(1:natoms) = zot
    atom_name(1:natoms) = atom_namet
    atom_type(1:natoms) = atom_typet
!    if (luselabels) then
!      catom_label(1:natoms) = catom_labelt
!    else
!      atom_label(1:natoms) = atom_labelt
!    end if
    atom_label(1:natoms) = atom_labelt
!    reference_number(1:natoms) = rnumt
!    reference_cell(1:natoms,:) = rcellt

    scale_factor = 3.07/2.35
    
    cell = cell*scale_factor
    a = a*scale_factor  ;  b = b*scale_factor  ;  c = c*scale_factor

    jatom = natoms
    do i = 1,natoms
      buffer = www_neighbours(i)
      read(buffer,*) j,neighbourlist
      neighbour_list: do j = 1,4
        if (i>neighbourlist(j)) cycle neighbour_list
        k = neighbourlist(j) + 1
        jatom = jatom + 1
        dx = xf(k) - xf(i) + 1.5d0
        dx = dx - aint(dx) - 0.5d0
        dy = yf(k) - yf(i) + 1.5d0
        dy = dy - aint(dy) - 0.5d0
        dz = zf(k) - zf(i) + 1.5d0
        dz = dz - aint(dz) - 0.5d0
        xf(jatom) = xf(i) + dx/2.0d0   
        yf(jatom) = yf(i) + dy/2.0d0   
        zf(jatom) = zf(i) + dz/2.0d0  
        xo(jatom) = xf(jatom)*a
        yo(jatom) = yf(jatom)*a
        zo(jatom) = zf(jatom)*a
        atom_name(jatom) = 'O '
        atom_type(jatom) = 8
      end do neighbour_list
    end do

    natoms = newatoms
    call assign_elements

    if (allocated(xft)) deallocate(xft,yft,zft)
    if (allocated(xot)) deallocate(xot,yot,zot)
    if (allocated(atom_namet)) deallocate(atom_namet)
    if (allocated(atom_typet)) deallocate(atom_typet)
    if (allocated(atom_labelt)) deallocate(atom_labelt)
    if (allocated(catom_labelt)) deallocate(catom_labelt)
    if (allocated(rnumt)) deallocate(rnumt)
    if (allocated(rcellt)) deallocate(rcellt)

    reference_number = 0
    reference_cell = 0


    if (ldiag) then
      write(main_output,*)
      write(main_output,*) 'Output from make_into_silica'
      write(main_output,*) ' Number of atoms: ',natoms
      write(main_output,*)
      write(main_output,'(a)') ' === Original configuration in fractional coordinates ==='
      do i = 1,natoms
        write(main_output,'(i0,2x,a4,2x,i3,3f10.4)') i,atom_name(i),atom_type(i),xf(i),yf(i),zf(i)
      end do
      write(main_output,*)
      write(main_output,'(a)') ' === Original configuration in orthogonal coordinates ==='
      do i = 1,natoms
        write(main_output,'(i0,2x,a4,2x,i3,3f10.4)') i,atom_name(i),atom_type(i),xo(i),yo(i),zo(i)
      end do
      write(main_output,*)
    end if



    return
  
  end subroutine make_into_silica


!===============================================================================


   subroutine calculate_pdf
!  ========================

!  Calculates the pair distribution function from either the starting structure or
!  the generated configuration.
!  There are several cases, and I will do each one of them in turn. These are
!  a) if we compute this using only a single unit cell
!  b) if we compute for a configuration
!  c) if that configuration is a nanoparticle

   use structure_data
   use arguments
   use fft_mod
   use channel_numbers

   implicit none

   double precision, parameter :: deltar = 0.02d0
   double precision, allocatable :: rpdf(:),pdfn(:),pdfx(:),dn(:),dx(:),florch(:),q(:), &
                                    qiqn(:),qiqx(:),sqn(:),sqx(:),t(:),w(:),partial(:,:,:), &
                                    g(:,:,:),gt(:)
   integer, parameter :: nsearch = 100
   integer :: nwidth,i,j,k,kk,npdfpts,ix,iy,iz,midpoint,n,ir,ik,jatom
   integer :: na,nb,nc,namax,nbmax,ncmax,one,nqiq,np,npairs,nlist
   integer, allocatable :: ip(:),ng(:,:,:),originlist(:)
   double precision, allocatable :: peak(:)
   double precision :: volume,afactor,pi,gxinfinity,gninfinity,deltaq,gscale
   double precision :: sigma, axb1,axb2,axb3,bxc1,bxc2,bxc3,cxa1,cxa2,cxa3,d,d1,d2,d3
   double precision :: rx,ry,rz,r,rr,x,gg
   double precision :: f,ftot,roi,roj,xomean,yomean,zomean,xjunk
   character(len=132) :: fname,buffer
   character(len=2), allocatable :: originatom(:)
   logical,allocatable :: luse(:)

   allocate(luse(natoms))
   luse = .true.

   if (llistfile) then
     luse = .false.
     open(ilist,file=trim(listfile),status='old',form='formatted')
       read(ilist,'(a)') buffer
       if (buffer(1:21)=='Number of particles =') then  ! We have reduced size AtomEye file
         read(buffer(22:),*) nlist
         do i = 1,10  ;  read(ilist,*)  ; end do
         allocate(originlist(nlist),originatom(nlist))
         do i = 1,nlist
           read(ilist,'(a)') buffer
           buffer = adjustl(buffer)
           n = index(buffer,' ')
           buffer = (adjustl(buffer(n:)))
           originatom(i) = buffer(1:2)
           buffer = adjustl(buffer(3:))
           read(buffer,*) xjunk,xjunk,xjunk,originlist(i)
           luse(originlist(i)) = .true.
         end do
       end if
     close(ilist)
   end if
   
   write(6,'(a,i0)') 'Number of atoms identified in origins file = ',count(luse)

   pi = acos(-1.0d0)
   if (lpdfbroaden) then
     sigma = pdfwidth
     nwidth = 4*pdfwidth/deltar
     allocate (peak(-nwidth:nwidth))
     do i = -nwidth,nwidth
       x = i*deltar/sigma
       peak(i) = exp(-x**2)/(sigma*sqrt(pi))
     end do
   end if   

   volume = cell(1,1)*cell(2,2)*cell(3,3) + cell(1,2)*cell(2,3)*cell(3,1) + &
            cell(1,3)*cell(2,1)*cell(3,2) - cell(1,1)*cell(2,3)*cell(3,2) - &
            cell(1,3)*cell(2,2)*cell(3,1) - cell(1,2)*cell(2,1)*cell(3,3)
   density = dble(natoms)/volume
   gscale = 4.0d0*pi*density*deltar

  if (lzeta) then
    xomean = sum(xo)/dble(natoms) ; yomean = sum(yo)/dble(natoms) ; zomean = sum(zo)/dble(natoms)
    xo = xo - xomean ; yo = yo - yomean ; zo = zo - zomean
  end if
   
!  Perform calculation for the nanoparticle   
   if (lnano) then
     if (abs(rpdfmax-max(rnano(1),rnano(2),rnano(3)))<deltar) then 
       rpdfmax = rpdfmax + deltar*nwidth
     else
       rpdfmax = min(rpdfmax,(max(rnano(1),rnano(2),rnano(3))+deltar*nwidth))
     end if
     npdfpts = nint(rpdfmax/deltar)
     allocate (rpdf(npdfpts),pdfn(npdfpts),pdfx(npdfpts))
     do i = 1,npdfpts
       rpdf(i) = dble(i)*deltar
     end do
     pdfn = 0.0d0
     pdfx = 0.0d0
     do i = 1,natoms
       if (lxyzlimits.and.(xf(i)<xyzlimits(1))) cycle
       if (lxyzlimits.and.(xf(i)>xyzlimits(2))) cycle
       if (lxyzlimits.and.(yf(i)<xyzlimits(3))) cycle
       if (lxyzlimits.and.(yf(i)>xyzlimits(4))) cycle
       if (lxyzlimits.and.(zf(i)<xyzlimits(5))) cycle
       if (lxyzlimits.and.(zf(i)>xyzlimits(6))) cycle
       do j = i+1,natoms
         rx = xo(i) - xo(j)
         ry = yo(i) - yo(j)
         rz = zo(i) - zo(j)
         r = sqrt(rx**2 + ry**2 + rz**2)
         if (r>rpdfmax) cycle
         midpoint = nint(r/deltar)
         peak = 0.0d0
         do k = -nwidth,nwidth
           kk = k + midpoint
           if (kk>npdfpts) cycle
           rr = rpdf(kk) - r
           peak(k) = exp(-(rr/sigma)**2)
           pdfn(kk) = pdfn(kk) + 2.0d0*peak(k)*blength(atom_type(i))*blength(atom_type(j))
           pdfx(kk) = pdfx(kk) + 2.0d0*peak(k)*abs(atom_type(i))*abs(atom_type(j))
         end do
       end do
     end do

!  Perform calculation based on one unit cell
   else if (lone.and.(.not.lconfigin).and.(.not.lonecell)) then
     write(6,'(a)') 'Calculating PDF for a single unit cell crystal structure'

   ! First check the grid we need using a random walk method
     namax = int(a/rpdfmax) + 1
     nbmax = int(b/rpdfmax) + 1
     ncmax = int(c/rpdfmax) + 1
     do n = 1,nsearch
       na = 0 ; nb = 0 ; nc = 0
       r = 0.0d0
       do while (rpdfmax>r)
         call random_number(x)
         i = int(x*3.0d0 + 1.0d0)
         call random_number(x)
         x = 2.0d0*x - 1.0d0
         one = 1 ; if (x<0) one = -1
         if (i==1) na = na + one
         if (i==2) nb = nb + one
         if (i==3) nc = nc + one
         rx = na*cell(1,1) + nb*cell(1,2) + nc*cell(1,3)
         ry = na*cell(2,1) + nb*cell(2,2) + nc*cell(2,3)
         rz = na*cell(3,1) + nb*cell(3,2) + nc*cell(3,3)
         r = sqrt(rx**2 + ry**2 + rz**2)
       end do
       if (abs(na)>namax) namax = abs(na)
       if (abs(nb)>nbmax) nbmax = abs(nb)
       if (abs(nc)>ncmax) ncmax = abs(nc)
     end do
   ! Now loop over the atoms
     npdfpts = nint(rpdfmax/deltar)
     allocate(rpdf(npdfpts),pdfn(npdfpts),pdfx(npdfpts))
     allocate(ng(npdfpts,ntypes,ntypes)) ; ng = 0
     allocate(g(npdfpts,ntypes,ntypes))
     do i = 1,npdfpts
       rpdf(i) = dble(i)*deltar
     end do
     pdfn = 0.0d0
     pdfx = 0.0d0
     do ix = -namax,namax
       do iy = -nbmax,nbmax
         do iz = -ncmax,ncmax
           do i = 1,natoms
             do j = 1,natoms
             if ((i==j).and.(ix==0).and.(iy==0).and.(iz==0)) cycle
               rx = xo(i) - xo(j) + ix*cell(1,1) + iy*cell(1,2) + iz*cell(1,3)
               ry = yo(i) - yo(j) + ix*cell(2,1) + iy*cell(2,2) + iz*cell(2,3)
               rz = zo(i) - zo(j) + ix*cell(3,1) + iy*cell(3,2) + iz*cell(3,3)
               r = sqrt(rx**2 + ry**2 + rz**2)
               if (r>rpdfmax) cycle
               midpoint = nint(r/deltar)
               kk = midpoint
               ng(kk,ordertype(i),ordertype(j)) = ng(kk,ordertype(i),ordertype(j)) + 1
!                if (lpdfbroaden) then
!                  peak = 0.0d0
!                  do k = -nwidth,nwidth
!                    kk = k + midpoint
!                    if (kk>npdfpts) cycle
!                    rr = rpdf(kk) - r
!                    peak(k) = exp(-(rr/sigma)**2)
!                    pdfn(kk) = pdfn(kk) + 2.0d0*peak(k)*blength(atom_type(i))*blength(atom_type(j))
!                    pdfx(kk) = pdfx(kk) + 2.0d0*peak(k)*abs(atom_type(i))*abs(atom_type(j))
!                  end do
!                else
!                  kk = midpoint
!                  pdfn(kk) = pdfn(kk) + 2.0d0*blength(atom_type(i))*blength(atom_type(j))
!                  pdfx(kk) = pdfx(kk) + 2.0d0*abs(atom_type(i))*abs(atom_type(j))
!                end if
             end do
           end do
         end do
       end do
     end do

!  Perform calculation based on one unit cell
   else if (lone.and.(.not.lconfigin).and.lonecell) then
     write(6,'(a)') 'Calculating PDF for a single unit cell only'

   ! Now Loop over the atoms
     npdfpts = nint(rpdfmax/deltar)
     allocate(rpdf(npdfpts),pdfn(npdfpts),pdfx(npdfpts))
     allocate(ng(npdfpts,ntypes,ntypes)) ; ng = 0
     allocate(g(npdfpts,ntypes,ntypes))
     do i = 1,npdfpts
       rpdf(i) = dble(i)*deltar
     end do
     pdfn = 0.0d0
     pdfx = 0.0d0
     do i = 1,natoms
       do j = 1,natoms
       if (i==j) cycle
         rx = xo(i) - xo(j)
         ry = yo(i) - yo(j)
         rz = zo(i) - zo(j)
         r = sqrt(rx**2 + ry**2 + rz**2)
         if (r>rpdfmax) cycle
         midpoint = nint(r/deltar)
         kk = midpoint
         ng(kk,ordertype(i),ordertype(j)) = ng(kk,ordertype(i),ordertype(j)) + 1
       end do
     end do

!  Perform calculation for the configuration
   else if (lconfigin) then
     write(6,'(a)') 'Calculating PDF for a configuration'
     volume = abs(cell(1,1)*cell(2,2)*cell(3,3) + cell(2,1)*cell(3,2)*cell(1,3) + &
                  cell(3,1)*cell(1,2)*cell(2,3) - cell(3,1)*cell(2,2)*cell(1,3) - &
                  cell(2,1)*cell(1,2)*cell(3,3) - cell(1,1)*cell(3,2)*cell(2,3))
     axb1 = cell(2,1)*cell(3,2) - cell(3,1)*cell(2,2)
     axb2 = cell(3,1)*cell(1,2) - cell(1,1)*cell(3,2)
     axb3 = cell(1,1)*cell(2,2) - cell(2,1)*cell(1,2)
     bxc1 = cell(2,2)*cell(3,3) - cell(3,2)*cell(2,3)
     bxc2 = cell(3,2)*cell(1,3) - cell(1,2)*cell(3,3)
     bxc3 = cell(1,2)*cell(2,3) - cell(2,2)*cell(1,3)
     cxa1 = cell(2,3)*cell(3,1) - cell(3,3)*cell(2,1)
     cxa2 = cell(3,3)*cell(1,1) - cell(1,3)*cell(3,1)
     cxa3 = cell(1,3)*cell(2,1) - cell(2,3)*cell(1,1)
     d1 = 1.0d0/dsqrt(axb1**2+axb2**2+axb3**2)
     d2 = 1.0d0/dsqrt(bxc1**2+bxc2**2+bxc3**2)
     d3 = 1.0d0/dsqrt(cxa1**2+cxa2**2+cxa3**2)
     d = volume/min(d1,d2,d3)
     rpdfmax = min(d,rpdfmax) 
     npdfpts = nint(rpdfmax/deltar)
!      if (ntypes==0) then
!        do i = -1,nelements
!          if (n_elements(i)>0) ntypes = ntypes + 1
!        end do
!      end if
!     if (.not.allocated(numoftype)) then
!       allocate(numoftype(ntypes))
!       numoftype = 0
!       i = 0
!       do j = -1,nelements
!         if (n_elements(j)>0) then
!           i = i + 1
!           numoftype(i) = n_elements(j)
!         end if
!       end do
!     end if

     npdfpts = nint(rpdfmax/deltar)
     allocate(rpdf(npdfpts),pdfn(npdfpts),pdfx(npdfpts))
     allocate(ng(npdfpts,ntypes,ntypes)) ; ng = 0
     allocate(g(npdfpts,ntypes,ntypes))
     do i = 1,npdfpts
       rpdf(i) = dble(i)*deltar
     end do

     pdfn = 0.0d0
     pdfx = 0.0d0
     do ix = -1,1
       do iy = -1,1
         do iz = -1,1
           do i = 1,natoms
            if (.not.luse(i)) cycle
            if (lxyzlimits.and.(xf(i)<xyzlimits(1))) cycle
            if (lxyzlimits.and.(xf(i)>xyzlimits(2))) cycle
            if (lxyzlimits.and.(yf(i)<xyzlimits(3))) cycle
            if (lxyzlimits.and.(yf(i)>xyzlimits(4))) cycle
            if (lxyzlimits.and.(zf(i)<xyzlimits(5))) cycle
            if (lxyzlimits.and.(zf(i)>xyzlimits(6))) cycle
            do jatom = 1,neighbours(i,0)
             j = neighbours(i,jatom)
             if ((i==j).and.(ix==0).and.(iy==0).and.(iz==0)) cycle
               rx = xo(i) - xo(j) + ix*cell(1,1) + iy*cell(1,2) + iz*cell(1,3)
               ry = yo(i) - yo(j) + ix*cell(2,1) + iy*cell(2,2) + iz*cell(2,3)
               rz = zo(i) - zo(j) + ix*cell(3,1) + iy*cell(3,2) + iz*cell(3,3)
               r = sqrt(rx**2 + ry**2 + rz**2)
               if (r>rpdfmax) cycle
               midpoint = nint(r/deltar)
               kk = midpoint
               ng(kk,ordertype(i),ordertype(j)) = ng(kk,ordertype(i),ordertype(j)) + 1
!                if (lpdfbroaden) then
!                  peak = 0.0d0
!                  do k = -nwidth,nwidth
!                    kk = k + midpoint
!                    if (kk>npdfpts) cycle
!                    rr = rpdf(kk) - r
!                    peak(k) = exp(-(rr/sigma)**2)
!                    pdfn(kk) = pdfn(kk) + 2.0d0*peak(k)*blength(atom_type(i))*blength(atom_type(j))
!                    pdfx(kk) = pdfx(kk) + 2.0d0*peak(k)*abs(atom_type(i))*abs(atom_type(j))
!                  end do
!                else
!                  kk = midpoint
!                  if (.not.lzeta) then
!                    pdfn(kk) = pdfn(kk) + 2.0d0*blength(atom_type(i))*blength(atom_type(j))
!                    pdfx(kk) = pdfx(kk) + 2.0d0*abs(atom_type(i))*abs(atom_type(j))
!                  else
!                    roi = sqrt(xo(i)**2 + yo(i)**2 + zo(i)**2)
!                    roj = sqrt(xo(j)**2 + yo(j)**2 + zo(j)**2)
!                    f = exp(-zeta*(roi+roj))
!                    ftot = ftot + f
!                    npairs = npairs + 1
!                    pdfn(kk) = pdfn(kk) + 2.0d0*f*blength(atom_type(i))*blength(atom_type(j))
!                    pdfx(kk) = pdfx(kk) + 2.0d0*f*abs(atom_type(i))*abs(atom_type(j))
!                  end if
!                end if
             end do
           end do
         end do
       end do
     end do
   else
     return
   end if

     do i = 1,ntypes
       g(:,i,:) = dble(ng(:,i,:))/dble(numoftype(i))
     end do
     do i = 1,npdfpts
       g(i,:,:) = g(i,:,:)/rpdf(i)**2
     end do
     g = g/gscale
     do j = 1,ntypes
       g(:,:,j) = g(:,:,j)/concentration(j)
     end do

  if (lpdfbroaden) then
    allocate(gt(npdfpts))
    do i = 1,ntypes
      do j = 1,ntypes
        gt = 0.0d0
        do ir = 1,npdfpts
          do k = -nwidth,nwidth
            ik = ir + k
            if (ik<1) then
              gg = 0.0d0
            else if (ik>npdfpts) then
              gg = 1.0d0
            else
              gg = g(ik,i,j)
            end if
            gt(ir) = gt(ir) + gg*peak(k)*deltar
          end do
        end do  
        g(:,i,j) = gt
      end do
    end do
  end if

  do i = 1,ntypes
    do j = 1,ntypes
      n = index(fileroot,'.')
      fname = trim(wfolder)//fileroot(1:n-1)//'_gr_'
      fname = trim(fname)//trim(element_of_type(i))//'_'//trim(element_of_type(j))//'.csv'
      open(ipdfn,file=trim(fname),form='formatted',status='unknown')
      do ir = 1,npdfpts
        write(ipdfn,'(f8.3,a,f14.3,a,f14.3)') rpdf(ir),' , ',g(ir,i,j),' , ',rpdf(ir)*(g(ir,i,j)-1.0d0)
      end do
      close(ipdfn)
    end do
  end do
  
  stop

   do i = 1,npdfpts
     pdfn(i) = pdfn(i)/rpdf(i)**2 
     pdfx(i) = pdfx(i)/rpdf(i)**2 
   end do
   if (lpdfbroaden) then
     afactor = 1.0d0/(sigma*sqrt(pi))
     afactor = afactor/(4.0d0*pi*density)/2.0d0 ! Factor of 2 needed to offset previous factor of 2
   else
     afactor = 1.0d0/(4.0d0*pi*density)
   end if
   pdfn = pdfn*afactor/dble(natoms)
   pdfx = pdfx*afactor/dble(natoms)

   n = index(fileroot,'.')
   if (n>0) then
     open(ipdfn,file=trim(wfolder)//fileroot(1:n-1)//'_pdfn.csv',form='formatted',status='unknown')
     open(ipdfx,file=trim(wfolder)//fileroot(1:n-1)//'_pdfx.csv',form='formatted',status='unknown')
   else
     open(ipdfn,file=trim(wfolder)//trim(fileroot)//'_pdfn.csv',form='formatted',status='unknown')
     open(ipdfx,file=trim(wfolder)//trim(fileroot)//'_pdfx.csv',form='formatted',status='unknown')
   end if

   do i = 1,npdfpts
     write(ipdfn,'(f8.3,a,f14.3)') rpdf(i),' , ',pdfn(i)
     write(ipdfx,'(f8.3,a,f14.3)') rpdf(i),' , ',pdfx(i)
   end do
   close(ipdfn)
   close(ipdfx)

   allocate(dx(npdfpts),dn(npdfpts),florch(npdfpts))
   gxinfinity = 0.0d0 ; gninfinity = 0.0d0
   do i = 1,natoms
     gxinfinity = gxinfinity + atom_type(i)
     gninfinity = gninfinity + blength(atom_type(i))
   end do
   gxinfinity = (gxinfinity/dble(natoms))**2 ; gninfinity = (gninfinity/dble(natoms))**2
   dx = pdfx - gxinfinity ; dx = dx*rpdf
   dn = pdfn - gninfinity ; dn = dn*rpdf
   do i = -nwidth,0
     dx(npdfpts+i) = 0.0d0
     dn(npdfpts+i) = 0.0d0
   end do

  if (lzeta) then
    pdfn = pdfn*dble(npairs)/ftot
    pdfx = pdfx*dble(npairs)/ftot
  end if

   n = index(fileroot,'.')
   if (n>0) then
     open(ipdfn,file=trim(wfolder)//fileroot(1:n-1)//'_pdfn.csv',form='formatted',status='unknown')
     open(ipdfx,file=trim(wfolder)//fileroot(1:n-1)//'_pdfx.csv',form='formatted',status='unknown')
   else
     open(ipdfn,file=trim(wfolder)//trim(fileroot)//'_pdfn.csv',form='formatted',status='unknown')
     open(ipdfx,file=trim(wfolder)//trim(fileroot)//'_pdfx.csv',form='formatted',status='unknown')
   end if

   write(ipdfn,'(a)') 'Distance (Ang) ,  Neutron G(r) ,  Neutron D(r)  ,  Neutron n(r)'
   write(ipdfx,'(a)') 'Distance (Ang) ,  X-ray G(r)   ,  X-ray D(r)    , X-ray n(r)'
   do i = 1,npdfpts
     write(ipdfn,'(f8.3,3(a,f14.4))') rpdf(i),' , ',pdfn(i),' , ',dn(i),' , ',sum(pdfn(1:i))
     write(ipdfx,'(f8.3,3(a,f14.4))') rpdf(i),' , ',pdfx(i),' , ',dx(i),' , ',sum(pdfx(1:i))
   end do
   close(ipdfn)
   close(ipdfx)


   do i = 1,npdfpts
     florch(i) = sin(pi*rpdf(i)/rpdfmax)/(pi*rpdf(i)/rpdfmax)
   end do
   dx = dx*florch ; dn = dn*florch
   np = int(log(dble(npdfpts))/log(2.0d0)) + 1
   nqiq = 2**np
   allocate(qiqn(0:nqiq-1),qiqx(0:nqiq-1),sqn(0:nqiq-1),sqx(0:nqiq-1),q(0:nqiq-1))
   qiqn = 0.0d0 ; qiqx = 0.0d0
   qiqn(1:npdfpts) = dn ; qiqx(1:npdfpts) = dx
   allocate(t(0:((nqiq/2)-1)),ip(0:int((2+sqrt(dble(nqiq/4))))),w(0:(nqiq*5/8-1)))
   ip = 0 ; t = 0.0d0 ; w = 0.0d0 
   call dfst(nqiq,qiqn,t,ip,w)
   ip = 0 ; t = 0.0d0 ; w = 0.0d0 
   call dfst(nqiq,qiqx,t,ip,w)   

   n = index(fileroot,'.')
   if (n>0) then
     open(ipdfn,file=trim(wfolder)//fileroot(1:n-1)//'_sqn.csv',form='formatted',status='unknown')
     open(ipdfx,file=trim(wfolder)//fileroot(1:n-1)//'_sqx.csv',form='formatted',status='unknown')
   else
     open(ipdfn,file=trim(wfolder)//trim(fileroot)//'_sqn.csv',form='formatted',status='unknown')
     open(ipdfx,file=trim(wfolder)//trim(fileroot)//'_sqx.csv',form='formatted',status='unknown')
   end if

   deltaq = 2.0d0*pi/dble(nqiq*deltar)
   do i = 0,nqiq-1
     q(i) = i*deltaq/2.0d0   ! Seemingly I need the factor of 2
   end do

   write(ipdfn,'(a)') 'Q (1/Ang) ,  Neutron Qi(Q) ,  Neutron i(Q)'
   write(ipdfx,'(a)') 'Q (1/Ang) ,  X-ray Qi(Q) ,  X-ray i(Q)'

   sqn(1:) = qiqn(1:)/q(1:) ; sqx(1:) = qiqx(1:)/q(1:)
   do i = 1,nqiq-1
     write(ipdfn,'(f8.3,2(a,f14.4))') q(i),' , ',qiqn(i),' , ',sqn(i)
     write(ipdfx,'(f8.3,2(a,g14.6))') q(i),' , ',qiqx(i),' , ',sqx(i)
   end do
   close(ipdfn)
   close(ipdfx)

   deallocate(rpdf,pdfn,pdfx,dx,dn,florch,qiqn,qiqx,q,sqn,sqx)

  if (lzeta) then
    xo = xo + xomean ; yo = yo + yomean ; zo = zo + zomean
  end if
   
   return
   
   end subroutine calculate_pdf
   

!===============================================================================
   
   
  subroutine dlpoly_pdf
! =====================
  
   use structure_data
   use arguments
   use fft_mod
   use channel_numbers

   implicit none

   integer :: nsets,npdfpts,i,j,k,n
   double precision, allocatable :: rpdf(:),pdfn(:),pdfx(:),g(:,:),factorn(:),factorx(:)
   double precision :: gscale,pi,conci,concj
   character(len=20), allocatable :: catompair(:)
   character(len=2), allocatable :: ci(:),cj(:)
   character(len=20) :: ctemp
   integer, allocatable :: nelementi(:),nelementj(:)
   
  pi = acos(-1.00)
  
  open(irdfdat,file='RDFDAT',form='formatted',status='old')
  
  read(irdfdat,*)
  read(irdfdat,*) nsets,npdfpts

  if (nsets/=(ntypes*(ntypes+1)/2)) then
    write(6,'(a)') 'There is a big inconsistency: the number of partial PDFs in RDFDAT does not'
    write(6,'(a)') 'match the number expected based on the number expected from the configuration'
    write(6,'(a)') 'file'
    write(6,'(a,i0)') 'The number of atom types = ',ntypes
    write(6,'(a,i0)') 'So number of partial PDFs expected = ',ntypes*(ntypes+1)/2
    write(6,'(a,i0)') 'But the number of partial PDFs in RDFDAT = ',nsets
    stop
  end if
  
  allocate(rpdf(npdfpts),pdfn(npdfpts),pdfx(npdfpts))
  allocate(g(npdfpts,nsets))
  allocate(catompair(nsets))
  allocate(nelementi(nsets),nelementj(nsets))
  allocate(factorn(nsets),factorx(nsets))
  allocate(ci(nsets),cj(nsets))

!    volume = cell(1,1)*cell(2,2)*cell(3,3) + cell(1,2)*cell(2,3)*cell(3,1) + &
!             cell(1,3)*cell(2,1)*cell(3,2) - cell(1,1)*cell(2,3)*cell(3,2) - &
!             cell(1,3)*cell(2,2)*cell(3,1) - cell(1,2)*cell(2,1)*cell(3,3)
!    density = dble(natoms)/volume
   gscale = 4.0d0*pi*density

  do n = 1,nsets
    read(irdfdat,'(a)') catompair(n)
    catompair(n) = adjustl(catompair(n))
    do i = 1,npdfpts
      read(irdfdat,*) rpdf(i),g(i,n)
    end do
  end do
  close(irdfdat)

  do n = 1,nsets
    ctemp = catompair(n)
    ci(n) = ctemp(1:2)
    ctemp = adjustl(ctemp(3:))
    cj(n) = ctemp(1:2)
  end do
  
  do n = 1,nsets
    do i = -1,nelements
      if (ci(n)==element(i)) then
        nelementi(n) = i
        cycle
      end if
    end do
  end do
  do n = 1,nsets
    do i = -1,nelements
      if (cj(n)==element(i)) then
        nelementj(n) = i
        cycle
      end if
    end do
  end do
  do n = 1,nsets
    i = nelementi(n) ; j = nelementj(n)
    conci = dble(n_elements(i))/dble(natoms) ; concj = dble(n_elements(j))/dble(natoms)
    factorn(n) = blength(i)*blength(j)*conci*concj
    factorx(n) = dble(i*j)*conci*concj
    if (i/=j) then
      factorn(n) = 2.0d0*factorn(n)
      factorx(n) = 2.0d0*factorx(n)
    end if
    write(6,*) i,j,n_elements(i),n_elements(j),conci,concj,factorn(n),factorx(n)
  end do
  
  g = g - 1.0d0  
  pdfn = 0.0d0 ; pdfx = 0.0d0
  
  do n = 1,nsets
    do i = 1,npdfpts
      pdfn(i) = pdfn(i) + factorn(n)*g(i,n)*rpdf(i)
      pdfx(i) = pdfx(i) + factorx(n)*g(i,n)*rpdf(i)
    end do
  end do

!   do i = 1,nsets
!     write(6,'(a)') 'Set 1 has atom types '//trim(ci(i))//' and '//trim(cj(i))//' with fractions and neutron/xray factors'
!     write(6,'(i2,2x,i2,2x,4f10.4)') nelementi(i),nelementj(i), &
!                          dble(n_elements(nelementi(i)))/dble(natoms), &
!                          dble(n_elements(nelementj(i)))/dble(natoms), &
!                         factorn(i),factorx(i)
!   end do
  
  pdfn = gscale*pdfn
  pdfx = gscale*pdfx

  open(ipdfn,file='DLPOLY_npdf.csv',form='formatted',status='unknown')
  open(ipdfx,file='DLPOLY_xpdf.csv',form='formatted',status='unknown')
  do i = 1,npdfpts
    write(ipdfn,'(f8.4,a,f12.6)') rpdf(i),' , ',pdfn(i)
    write(ipdfx,'(f8.4,a,f12.6)') rpdf(i),' , ',pdfx(i)
  end do
  close(ipdfn)
  close(ipdfx)
  

  return
  
  end subroutine dlpoly_pdf
  
  !===============================================================================

  subroutine assign_element_Uiso
  
  use structure_data
  use arguments
  
  implicit none
  
  integer :: n, i, ierror, itype, nUiso
  character(len=2) :: ctype 
  character(len=300) :: corder_list_temp, cUiso_list_temp, catomU
  double precision :: cUiso
  
  Uiso = 0.0 
  nUiso = 0
  corder_list_temp = corder_list
  cUiso_list_temp = cUiso_list
  
  if (.not. lorder_list) then
      write(6,'(a)') '-Uiso or -Biso option can only be used with -order option'
      stop
  endif
  
  call remove_spaces(corder_list_temp)
  
  do while (len_trim(corder_list_temp) > 0)
      n = index(corder_list_temp, ' ')
      if (n > 0) then
          corder_list_temp = corder_list_temp(n+1:)
          nUiso = nUiso + 1
      endif
  enddo
  
  corder_list_temp = corder_list
  
  
  if (.not. lUiso_list) then
      itype = 0
      cUiso_list = ''
      write(6,'(a)') 'Please give Uiso for each atom type: '
      atom_list_loop_U: do while (itype<nUiso)
        itype = itype + 1
        write(6,'(a,i0,a)', advance="no") 'Uiso of atom for order position ',itype,': '
        read(5,'(a)') catomU
        cUiso_list = trim(cUiso_list) // ' ' // trim(catomU)
      enddo atom_list_loop_U
  endif
  
  cUiso_list_temp = cUiso_list
  
  ierror = 0
  read(cUiso_list_temp, *, iostat=ierror) cUiso
  n = index(cUiso_list_temp, ' ')
  cUiso_list_temp = cUiso_list_temp(n+1:)
  n = index(corder_list_temp, ' ')
  ctype = ''
  ctype(1:n-1) = corder_list_temp(1:n-1)
  corder_list_temp = corder_list_temp(n+1:)
  
  do_elements_loop_1: do i = -1, nelements
		  if ((ichar(ctype(1:1))>=ichar('a')).and.(ichar(ctype(1:1))<=ichar('z'))) &
                ctype(1:1) = char(ichar(ctype(1:1)) - ichar('a') + ichar('A'))
		  if ((ichar(ctype(2:2))>=ichar('A')).and.(ichar(ctype(2:2))<=ichar('Z'))) &
                ctype(2:2) = char(ichar(ctype(2:2)) - ichar('A') + ichar('a'))		
          if (element(i) == ctype) then
              Uiso(i) = cUiso
              exit do_elements_loop_1
          endif          
      end do do_elements_loop_1
  
  
  read(cUiso_list_temp, *, iostat=ierror) cUiso
  
  readloop1: do while (ierror/=-1)
      do_elements_loop: do i = -1, nelements
		  if ((ichar(ctype(1:1))>=ichar('a')).and.(ichar(ctype(1:1))<=ichar('z'))) &
                ctype(1:1) = char(ichar(ctype(1:1)) - ichar('a') + ichar('A'))
		  if ((ichar(ctype(2:2))>=ichar('A')).and.(ichar(ctype(2:2))<=ichar('Z'))) &
                ctype(2:2) = char(ichar(ctype(2:2)) - ichar('A') + ichar('a'))		
          if (element(i) == ctype) then
              Uiso(i) = cUiso
              exit do_elements_loop
          endif          
      end do do_elements_loop
      read(cUiso_list_temp, *, iostat=ierror) cUiso
      n = index(cUiso_list_temp, ' ')
      cUiso_list_temp = cUiso_list_temp(n+1:)
      n = index(corder_list_temp, ' ')
      ctype = ''
      ctype(1:n-1) = corder_list_temp(1:n-1)
      corder_list_temp = corder_list_temp(n+1:)
  enddo readloop1
  
  
  end subroutine assign_element_Uiso


!===============================================================================

  subroutine assign_element_names

! ===============================

  use structure_data

  integer :: i,j,lc,uc
  character(len=2) :: celement

    element(-1) = 'D ' ; element_name(-1) = 'Deuterium' ; element_mass(-1) = 2.01410178d0 ; blength(-1) = 0.66710d0
    element(0) = 'Va' ; element_name(0) = 'Vacancy' ;  element_mass(0) = 0.0d0 ; blength(0) = 0.0d0
    element(1) = 'H ' ; element_name(1) = 'Hydrogen' ;  element_mass(1) = 1.0794d0 ; blength(1) = -0.37390d0
    element(2) = 'He' ; element_name(2) = 'Helium' ;  element_mass(2) = 4.02602d0 ; blength(2) = 0.3260d0
    element(3) = 'Li' ; element_name(3) = 'Lithium' ;  element_mass(3) = 6.941d0 ; blength(3) = -0.190d0
    element(4) = 'Be' ; element_name(4) = 'Beryllium' ;  element_mass(4) = 9.012182d0 ; blength(4) = 0.7790d0
    element(5) = 'B ' ; element_name(5) = 'Boron' ;  element_mass(5) = 10.811d0 ; blength(5) = 0.530d0
    element(6) = 'C ' ; element_name(6) = 'Carbon' ;  element_mass(6) = 12.011d0 ; blength(6) = 0.66460d0
    element(7) = 'N ' ; element_name(7) = 'Nitrogen' ;  element_mass(7) = 14.0674d0 ; blength(7) = 0.9360d0
    element(8) = 'O ' ; element_name(8) = 'Oxygen' ;  element_mass(8) = 15.9994d0 ; blength(8) = 0.58030d0
    element(9) = 'F ' ; element_name(9) = 'Fluorine' ;  element_mass(9) = 18.9984032d0 ; blength(9) = 0.56540d0
    element(10) = 'Ne' ; element_name(10) = 'Neon' ;  element_mass(10) = 20.1797d0 ; blength(10) = 0.45660d0
    element(11) = 'Na' ; element_name(11) = 'Sodium' ;  element_mass(11) = 22.989768d0 ; blength(11) = 0.3580d0
    element(12) = 'Mg' ; element_name(12) = 'Magnesium' ;  element_mass(12) = 24.3050d0 ; blength(12) = 0.53750d0
    element(13) = 'Al' ; element_name(13) = 'Aluminum' ;  element_mass(13) = 26.981539d0 ; blength(13) = 0.34490d0
    element(14) = 'Si' ; element_name(14) = 'Silicon' ;  element_mass(14) = 28.0855d0 ; blength(14) = 0.415340d0
    element(15) = 'P ' ; element_name(15) = 'Phosphorus' ;  element_mass(15) = 30.973762d0 ; blength(15) = 0.5130d0
    element(16) = 'S ' ; element_name(16) = 'Sulfur' ;  element_mass(16) = 32.066d0 ; blength(16) = 0.28470d0
    element(17) = 'Cl' ; element_name(17) = 'Chlorine' ;  element_mass(17) = 35.4527d0 ; blength(17) = 0.95770d0
    element(18) = 'Ar' ; element_name(18) = 'Argon' ;  element_mass(18) = 39.948d0 ; blength(18) = 0.19090d0
    element(19) = 'K ' ; element_name(19) = 'Potassium' ;  element_mass(19) = 39.0983d0 ; blength(19) = 0.3670d0
    element(20) = 'Ca' ; element_name(20) = 'Calcium' ;  element_mass(20) = 40.078d0 ; blength(20) = 0.4760d0
    element(21) = 'Sc' ; element_name(21) = 'Scandium' ;  element_mass(21) = 44.955910d0 ; blength(21) = 1.2290d0
    element(22) = 'Ti' ; element_name(22) = 'Titanium' ;  element_mass(22) = 47.88d0 ; blength(22) = -0.34380d0
    element(23) = 'V ' ; element_name(23) = 'Vanadium' ;  element_mass(23) = 50.9415d0 ; blength(23) = -0.038240d0
    element(24) = 'Cr' ; element_name(24) = 'Chromium' ;  element_mass(24) = 51.9961d0 ; blength(24) = 0.36350d0
    element(25) = 'Mn' ; element_name(25) = 'Manganese' ;  element_mass(25) = 54.93805d0 ; blength(25) = -0.3730d0
    element(26) = 'Fe' ; element_name(26) = 'Iron' ;  element_mass(26) = 55.847d0 ; blength(26) = 0.9540d0
    element(27) = 'Co' ; element_name(27) = 'Cobalt' ;  element_mass(27) = 58.93320d0 ; blength(27) = 0.2780d0
    element(28) = 'Ni' ; element_name(28) = 'Nickel' ;  element_mass(28) = 58.6934d0 ; blength(28) = 1.030d0
    element(29) = 'Cu' ; element_name(29) = 'Copper' ;  element_mass(29) = 63.546d0 ; blength(29) = 0.77180d0
    element(30) = 'Zn' ; element_name(30) = 'Zinc' ;  element_mass(30) = 65.39d0 ; blength(30) = 0.5680d0
    element(31) = 'Ga' ; element_name(31) = 'Gallium' ;  element_mass(31) = 69.723d0 ; blength(31) = 0.72880d0
    element(32) = 'Ge' ; element_name(32) = 'Germanium' ;  element_mass(32) = 72.61d0 ; blength(32) = 0.81850d0
    element(33) = 'As' ; element_name(33) = 'Arsenic' ;  element_mass(33) = 74.92159d0 ; blength(33) = 0.6580d0
    element(34) = 'Se' ; element_name(34) = 'Selenium' ;  element_mass(34) = 78.96d0 ; blength(34) = 0.7970d0
    element(35) = 'Br' ; element_name(35) = 'Bromine' ;  element_mass(35) = 79.904d0 ; blength(35) = 0.67950d0
    element(36) = 'Kr' ; element_name(36) = 'Krypton' ;  element_mass(36) = 83.80d0 ; blength(36) = 0.7810d0
    element(37) = 'Rb' ; element_name(37) = 'Rubidium' ;  element_mass(37) = 85.4678d0 ; blength(37) = 0.7090d0
    element(38) = 'Sr' ; element_name(38) = 'Strontium' ;  element_mass(38) = 87.62d0 ; blength(38) = 0.7020d0
    element(39) = 'Y ' ; element_name(39) = 'Yttrium' ;  element_mass(39) = 88.90585d0 ; blength(39) = 0.7750d0
    element(40) = 'Zr' ; element_name(40) = 'Zirconium' ;  element_mass(40) = 91.224d0 ; blength(40) = 0.7160d0
    element(41) = 'Nb' ; element_name(41) = 'Niobium' ;  element_mass(41) = 92.90638d0 ; blength(41) = 0.70540d0
    element(42) = 'Mo' ; element_name(42) = 'Molybdenum' ;  element_mass(42) = 95.94d0 ; blength(42) = 0.67150d0
    element(43) = 'Tc' ; element_name(43) = 'Technetium' ;  element_mass(43) = 98.0d0 ; blength(43) = 0.680d0
    element(44) = 'Ru' ; element_name(44) = 'Ruthenium' ;  element_mass(44) = 101.07d0 ; blength(44) = 0.7210d0
    element(45) = 'Rh' ; element_name(45) = 'Rhodium' ;  element_mass(45) = 102.90550d0 ; blength(45) = 0.5880d0
    element(46) = 'Pd' ; element_name(46) = 'Palladium' ;  element_mass(46) = 106.42d0 ; blength(46) = 0.5910d0
    element(47) = 'Ag' ; element_name(47) = 'Silver' ;  element_mass(47) = 107.8682d0 ; blength(47) = 0.59220d0
    element(48) = 'Cd' ; element_name(48) = 'Cadmium' ;  element_mass(48) = 112.411d0 ; blength(48) = 0.510d0
    element(49) = 'In' ; element_name(49) = 'Indium' ;  element_mass(49) = 114.82d0 ; blength(49) = 0.40650d0
    element(50) = 'Sn' ; element_name(50) = 'Tin' ;  element_mass(50) = 118.710d0 ; blength(50) = 0.62250d0
    element(51) = 'Sb' ; element_name(51) = 'Antimony' ;  element_mass(51) = 121.757d0 ; blength(51) = 0.5570d0
    element(52) = 'Te' ; element_name(52) = 'Tellurium' ;  element_mass(52) = 127.60d0 ; blength(52) = 0.580d0
    element(53) = 'I ' ; element_name(53) = 'Iodine' ;  element_mass(53) = 126.90447d0 ; blength(53) = 0.5280d0
    element(54) = 'Xe' ; element_name(54) = 'Xenon' ;  element_mass(54) = 131.29d0 ; blength(54) = 0.4920d0
    element(55) = 'Cs' ; element_name(55) = 'Cesium' ;  element_mass(55) = 132.90543d0 ; blength(55) = 0.5420d0
    element(56) = 'Ba' ; element_name(56) = 'Barium' ;  element_mass(56) = 137.327d0 ; blength(56) = 0.5070d0
    element(57) = 'La' ; element_name(57) = 'Lanthanum' ;  element_mass(57) = 138.9055d0 ; blength(57) = 0.8240d0
    element(58) = 'Ce' ; element_name(58) = 'Cerium' ;  element_mass(58) = 140.115d0 ; blength(58) = 0.4840d0
    element(59) = 'Pr' ; element_name(59) = 'Praseodymium' ;  element_mass(59) = 140.90765d0 ; blength(59) = 0.4450d0
    element(60) = 'Nd' ; element_name(60) = 'Neodymium' ;  element_mass(60) = 144.24d0 ; blength(60) = 0.7690d0
    element(61) = 'Pm' ; element_name(61) = 'Promethium' ;  element_mass(61) = 145.0d0 ; blength(61) = 1.260d0
    element(62) = 'Sm' ; element_name(62) = 'Samarium' ;  element_mass(62) = 150.36d0 ; blength(62) = 0.080d0
    element(63) = 'Eu' ; element_name(63) = 'Europium' ;  element_mass(63) = 151.965d0 ; blength(63) = 0.7220d0
    element(64) = 'Gd' ; element_name(64) = 'Gadolinium' ;  element_mass(64) = 157.25d0 ; blength(64) = 0.650d0
    element(65) = 'Tb' ; element_name(65) = 'Terbium' ;  element_mass(65) = 158.92534d0 ; blength(65) = 0.7380d0
    element(66) = 'Dy' ; element_name(66) = 'Dysprosium' ;  element_mass(66) = 162.50d0 ; blength(66) = 1.690d0
    element(67) = 'Ho' ; element_name(67) = 'Holmium' ;  element_mass(67) = 164.93032d0 ; blength(67) = 0.8010d0
    element(68) = 'Er' ; element_name(68) = 'Erbium' ;  element_mass(68) = 167.26d0 ; blength(68) = 0.8160d0
    element(69) = 'Tm' ; element_name(69) = 'Thulium' ;  element_mass(69) = 168.93421d0 ; blength(69) = 0.7070d0
    element(70) = 'Yb' ; element_name(70) = 'Ytterbium' ;  element_mass(70) = 173.04d0 ; blength(70) = 1.2430d0
    element(71) = 'Lu' ; element_name(71) = 'Lutetium' ;  element_mass(71) = 174.967d0 ; blength(71) = 0.7210d0
    element(72) = 'Hf' ; element_name(72) = 'Hafnium' ;  element_mass(72) = 178.49d0 ; blength(72) = 0.7770d0
    element(73) = 'Ta' ; element_name(73) = 'Tantalum' ;  element_mass(73) = 180.9479d0 ; blength(73) = 0.6910d0
    element(74) = 'W ' ; element_name(74) = 'Tungsten' ;  element_mass(74) = 183.85d0 ; blength(74) = 0.4860d0
    element(75) = 'Re' ; element_name(75) = 'Rhenium' ;  element_mass(75) = 186.207d0 ; blength(75) = 0.920d0
    element(76) = 'Os' ; element_name(76) = 'Osmium' ;  element_mass(76) = 190.2d0 ; blength(76) = 1.070d0
    element(77) = 'Ir' ; element_name(77) = 'Iridium' ;  element_mass(77) = 192.22d0 ; blength(77) = 1.060d0
    element(78) = 'Pt' ; element_name(78) = 'Platinum' ;  element_mass(78) = 195.08d0 ; blength(78) = 0.960d0
    element(79) = 'Au' ; element_name(79) = 'Gold' ;  element_mass(79) = 196.96654d0 ; blength(79) = 0.7630d0
    element(80) = 'Hg' ; element_name(80) = 'Mercury' ;  element_mass(80) = 20.59d0 ; blength(80) = 1.26920d0
    element(81) = 'Tl' ; element_name(81) = 'Thallium' ;  element_mass(81) = 204.3833d0 ; blength(81) = 0.87760d0
    element(82) = 'Pb' ; element_name(82) = 'Lead' ;  element_mass(82) = 207.2d0 ; blength(82) = 0.94050d0
    element(83) = 'Bi' ; element_name(83) = 'Bismuth' ;  element_mass(83) = 208.98037d0 ; blength(83) = 0.85320d0
    element(84) = 'Po' ; element_name(84) = 'Polonium' ;  element_mass(84) = 209.0d0 ; blength(84) = 0.0d0
    element(85) = 'At' ; element_name(85) = 'Astatine' ;  element_mass(85) = 210.0d0 ; blength(85) = 0.0d0
    element(86) = 'Rn' ; element_name(86) = 'Radon' ;  element_mass(86) = 222.0d0 ; blength(86) = 0.0d0
    element(87) = 'Fr' ; element_name(87) = 'Francium' ;  element_mass(87) = 223.0d0 ; blength(87) = 0.0d0
    element(88) = 'Ra' ; element_name(88) = 'Radium' ;  element_mass(88) = 226.025d0 ; blength(88) = 1.0d0
    element(89) = 'Ac' ; element_name(89) = 'Actinium' ;  element_mass(89) = 227.028d0 ; blength(89) = 0.0d0
    element(90) = 'Th' ; element_name(90) = 'Thorium' ;  element_mass(90) = 232.0381d0 ; blength(90) = 1.0520d0
    element(91) = 'Pa' ; element_name(91) = 'Protactinium' ;  element_mass(91) = 231.0359d0 ; blength(91) = 0.910d0
    element(92) = 'U ' ; element_name(92) = 'Uranium' ;  element_mass(92) = 238.0289d0 ; blength(92) = 0.84170d0
    element(93) = 'Np' ; element_name(93) = 'Neptunium' ;  element_mass(93) = 237.048d0 ; blength(93) = 1.0550d0
    element(94) = 'Pu' ; element_name(94) = 'Plutonium' ;  element_mass(94) = 244.0d0 ; blength(94) = 0.0d0
    element(95) = 'Am' ; element_name(95) = 'Americium' ;  element_mass(95) = 243.0d0 ; blength(95) = 0.830d0
    element(96) = 'Cm' ; element_name(96) = 'Curium' ;  element_mass(96) = 247.0d0 ; blength(96) = 0.0d0
    element(97) = 'Bk' ; element_name(97) = 'Berkelium' ;  element_mass(97) = 247.0d0 ; blength(97) = 0.0d0
    element(98) = 'Cf' ; element_name(98) = 'Californium' ;  element_mass(98) = 251.0d0 ; blength(98) = 0.0d0
    element(99) = 'Es' ; element_name(99) = 'Einsteinium' ;  element_mass(99) = 252.0d0 ; blength(99) = 0.0d0
    element(10) = 'Fm' ; element_name(10) = 'Fermium' ;  element_mass(10) = 257.0d0 ; blength(10) = 0.0d0
    element(101) = 'Md' ; element_name(101) = 'Mendelevium' ;  element_mass(101) = 258.0d0 ; blength(101) = 0.0d0
    element(102) = 'No' ; element_name(102) = 'Nobelium' ;  element_mass(102) = 259.0d0 ; blength(102) = 0.0d0
    element(103) = 'Lr' ; element_name(103) = 'Lawrencium' ;  element_mass(103) = 262.0d0 ; blength(103) = 0.0d0
    element(104) = 'Rf' ; element_name(104) = 'Rutherfordium' ;  element_mass(104) = 267.0d0 ; blength(104) = 0.0d0
    element(105) = 'Db' ; element_name(105) = 'Dubnium' ;  element_mass(105) = 268.0d0 ; blength(105) = 0.0d0
    element(106) = 'Sg' ; element_name(106) = 'Seaborgium' ;  element_mass(106) = 271.0d0 ; blength(106) = 0.0d0
    element(107) = 'Bh' ; element_name(107) = 'Bohrium' ;  element_mass(107) = 270.0d0 ; blength(107) = 0.0d0
    element(108) = 'Hs' ; element_name(108) = 'Hassium' ;  element_mass(108) = 277.0d0 ; blength(108) = 0.0d0
    element(109) = 'Mt' ; element_name(109) = 'Meitnerium' ;  element_mass(109) = 276.0d0 ; blength(109) = 0.0d0
    element(110) = 'Ds' ; element_name(110) = 'Darmstadtium' ;  element_mass(110) = 281.0d0 ; blength(110) = 0.0d0

    rgb(1,:)  = (/ 1.00000 , 1.00000 , 1.00000 /)
    rgb(2,:)  = (/ 1.00000 , 0.78430 , 0.78430 /)
    rgb(3,:)  = (/ 0.64709 , 0.16469 , 0.16469 /)
    rgb(4,:)  = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(5,:)  = (/ 0.00000 , 1.00000 , 0.00000 /)
    rgb(6,:)  = (/ 0.60000 , 0.60000 , 0.60000 /)
    rgb(7,:)  = (/ 0.56080 , 0.56080 , 1.00000 /)
    rgb(8,:)  = (/ 0.94119 , 0.00000 , 0.00000 /)
    rgb(9,:)  = (/ 0.78430 , 0.64709 , 0.09409 /)
    rgb(10,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(11,:) = (/ 0.00000 , 0.00000 , 1.00000 /)
    rgb(12,:) = (/ 0.16469 , 0.50199 , 0.16469 /)
    rgb(13,:) = (/ 0.50199 , 0.50199 , 0.56469 /)
    rgb(14,:) = (/ 0.78430 , 0.64709 , 0.09409 /)
    rgb(15,:) = (/ 1.00000 , 0.64709 , 0.00000 /)
    rgb(16,:) = (/ 1.00000 , 0.78430 , 0.19609 /)
    rgb(17,:) = (/ 0.00000 , 1.00000 , 0.00000 /)
    rgb(18,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(19,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(20,:) = (/ 0.50199 , 0.50199 , 0.56469 /)
    rgb(21,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(22,:) = (/ 0.50199 , 0.50199 , 0.56469 /)
    rgb(23,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(24,:) = (/ 0.50199 , 0.50199 , 0.56469 /)
    rgb(25,:) = (/ 0.50199 , 0.50199 , 0.56469 /)
    rgb(26,:) = (/ 1.00000 , 0.64709 , 0.00000 /)
    rgb(27,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(28,:) = (/ 0.64709 , 0.16469 , 0.16469 /)
    rgb(29,:) = (/ 0.64709 , 0.16469 , 0.16469 /)
    rgb(30,:) = (/ 0.64709 , 0.16469 , 0.16469 /)
    rgb(31,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(32,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(33,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(34,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(35,:) = (/ 0.64709 , 0.16469 , 0.16469 /)
    rgb(36,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(37,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(38,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(39,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(40,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(41,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(42,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(43,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(44,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(45,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(46,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(47,:) = (/ 0.50199 , 0.50199 , 0.56469 /)
    rgb(48,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(49,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(50,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(51,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(52,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(53,:) = (/ 0.56080 , 0.56080 , 1.00000 /)
    rgb(54,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(55,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(56,:) = (/ 1.00000 , 0.64709 , 0.00000 /)
    rgb(57,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(58,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(59,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(60,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(61,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(62,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(63,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(64,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(65,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(66,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(67,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(68,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(69,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(70,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(71,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(72,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(73,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(74,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(75,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(76,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(77,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(78,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(79,:) = (/ 0.78430 , 0.64709 , 0.09409 /)
    rgb(80,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(81,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(82,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(83,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(84,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(85,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(86,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(87,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(88,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(89,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(90,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(91,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(92,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(93,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(94,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(95,:) = (/ 1.00000 , 0.07839 , 0.57649 /)
    rgb(96,:) = (/ 0.00000 , 0.00000 , 0.00000 /)
    rgb(97,:) = (/ 0.00000 , 0.00000 , 0.00000 /)
    rgb(98,:) = (/ 0.00000 , 0.00000 , 0.00000 /)
    rgb(99,:) = (/ 0.00000 , 0.00000 , 0.00000 /)
    rgb(10,:) = (/ 0.00000 , 0.00000 , 0.00000 /)
    rgb(101,:) = (/0.00000 , 0.00000 , 0.00000 /)
    rgb(102,:) = (/0.00000 , 0.00000 , 0.00000 /)
    rgb(103,:) = (/0.00000 , 0.00000 , 0.00000 /)
    rgb(104,:) = (/0.00000 , 0.00000 , 0.00000 /)
    rgb(105,:) = (/0.00000 , 0.00000 , 0.00000 /)
    rgb(106,:) = (/0.00000 , 0.00000 , 0.00000 /)
    rgb(107,:) = (/0.00000 , 0.00000 , 0.00000 /)
    rgb(108,:) = (/0.00000 , 0.00000 , 0.00000 /)
    rgb(109,:) = (/0.00000 , 0.00000 , 0.00000 /)
    rgb(110,:) = (/0.00000 , 0.00000 , 0.00000 /)

    do i = -1,nelements
      celement = element(i)
      if(len(celement)==2) then
        uc = ichar('A') ; lc = ichar('a')
        j = ichar(celement(2:2)) - lc + uc
        celement(2:2) = char(j)
        elementuc(i) = celement
      end if
    end do

    do i = -1,nelements
      celement = element(i)
      celement(1:1) = char( ichar(celement(1:1)) - ichar('A') + ichar('a') )
      elementlc(i) = celement
    end do

    return

  end subroutine assign_element_names

!===============================================================================

  subroutine assign_elements
! ==========================

    use structure_data

    implicit none

    integer :: i,j,iel
    logical :: lassigned
    character(len=4) :: cjunk

    n_elements = 0

    atom_loop: do i = 1,natoms
    lassigned = .false.
    element_loop: do while (.not.lassigned)
! ..... Loop for elements with two letters in name
      sub_loop1: do iel = -1,nelements
        if (len(trim(element(iel)))/=2) cycle sub_loop1
        if ((index(atom_name(i),trim(element(iel)))==1).or. &
            (index(atom_name(i),trim(elementuc(iel)))==1).or. &
            (index(atom_name(i),trim(elementlc(iel)))==1) ) then
          atom_type(i) = iel
          lassigned = .true.
          exit element_loop
        end if
      end do sub_loop1
! ..... Loop for elements with one letter in name
      sub_loop2: do iel = -1,nelements
        if (len(trim(element(iel)))/=1) cycle sub_loop2
        if ((index(atom_name(i),trim(element(iel)))==1).or. &
            (index(atom_name(i),trim(elementuc(iel)))==1).or. &
            (index(atom_name(i),trim(elementlc(iel)))==1) ) then
          atom_type(i) = iel
          lassigned = .true.
          exit element_loop
        end if
      end do sub_loop2
! ..... Ask if element can't be assigned
      cjunk = atom_name(i)
      write(6,'(3a)') 'Element type for atom label ',trim(atom_name(i)),' not recognised.'
      write(6,'(a,$)') 'Please give element symbol: '
      read(5,'(a)') atom_name(i)
        do j = i+1,natoms
          if (atom_name(j)==cjunk) atom_name(j) = atom_name(i)
        end do
        do j = 1,ntypes
          if (aname(j)==cjunk) aname(j) = atom_name(i)
        end do
    end do element_loop
    n_elements(atom_type(i)) = n_elements(atom_type(i)) + 1
    end do atom_loop

  end subroutine assign_elements

!===============================================================================

  subroutine read_error(message,file)
! ===================================

    implicit none

    character*(*) :: message,file

    write(6,'(a)') 'File read error in file of type '//file
    write(6,'(a)') 'At point of reading: '//message
    write(6,'(a)') 'Program ending; you need to look at the faulty file'

    stop

  end subroutine read_error

!===============================================================================

   subroutine order_error()
! ===================================

    implicit none

    write(6,'(a)') 'Program ending; you need to look at the -order[] parameter'
    write(6,'(a)') 'Remember also to include Va if necessary'

    stop

   end subroutine order_error

!===============================================================================
   
      subroutine symmetry_error()
! ===================================

    implicit none

    write(6,'(a)') 'Program ending. There is no symmetry operators read from the file'
    write(6,'(a)') 'In case of cif file a list of symmetry operators have to be provided'
    write(6,'(a)') 'The program does not read spacegroup name'

    stop

  end subroutine symmetry_error

!===============================================================================

  subroutine end_of_file_error(message,file)
! ==========================================

    implicit none

    character*(*) :: message,file

    write(6,'(a)') 'End-of-file error in file of type '//file
    write(6,'(a)') 'At point of reading: '//message
    write(6,'(a)') 'Program ending; you need to look at the faulty file'

    stop

  end subroutine end_of_file_error

!===============================================================================

  subroutine get_expfile_parameters
! =================================

!--------------------------------------------------------------------------
!
!     Subroutine to generate the instrument parameter file, background
!     parameter file, and the Bragg scattering data file for RMCProfile.
!
!     This is called if the structure file type is .TBL and the .exp
!     file exists.
!
!--------------------------------------------------------------------------

  use arguments
  use structure_data
  use annotations
  use channel_numbers

  implicit none

  integer :: ihist,nhist,ierror,i,j,ncheck, &
           junk,jtype,npts,n,nptsused,nhist2,n1,n2,ihist2,ihtyp,ipol
  integer, allocatable :: nback(:),bank(:),nbraggpts(:),ipt(:),histpts(:),bank2(:),hist2(:)
  double precision :: zzz,peak(12),difc,difa,zero,theta,scale,volume,tmin,tmax,lambda, polar,abs_corr
  double precision :: pi
  double precision, allocatable :: back(:),tof(:),ydata(:),tofh(:),ydatah(:)
  logical :: l1,l2,lst_exists
  logical, allocatable :: usehist(:)
  character(len=120) :: fileexp,fileinst,fileback,filelst,filebragg,dataname
  character(len=132), allocatable :: cprofile(:,:),cicons(:),chscale(:),cbnkpar(:),ctrnge(:),cabspar(:)
  character(len=1000), allocatable :: cback(:,:)
  character(len=132) :: buffer,chap,cbk,cpts, cab, buff_temp
  character(len=10) :: shortbuffer
  character(len=3), allocatable :: hist_type(:)
  integer :: abs_type

  pi = acos(-1.0d0)

  if (ldiag) then
    write(main_output,'(/a)') 'Output from get_expfile_parameters'
    write(main_output,'(a/)') '=================================='
  end if


  n = index(fileroot,'.')
  fileexp = filename(1:n)//'EXP'
  inquire(file=trim(fileexp),exist=l1)
  fileexp = filename(1:n)//'exp'
  inquire(file=trim(fileexp),exist=l2)
  if (l1) fileexp = trim(filename(1:n)//'EXP')
  if (l2) fileexp = trim(filename(1:n)//'exp')
  fileinst = trim(wfolder)//trim(fileroot(1:n)//'inst')

  open(iexp,file=trim(fileexp),form='formatted',status='old')

! Extract number of histograms
  histcount: do
    read(iexp,'(a)') buffer
    if (index(buffer,'EXPR  NHST')>0) then
      read(buffer(13:),*) nhist                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
      exit histcount
    end if
  end do histcount
  rewind(iexp)

  allocate(usehist(nhist))
  allocate(hist_type(nhist))
  if(allocated(useonly)) then ! If we use the -usehist option, only process the marked histograms
     usehist = .false.
     usehist(useonly) = .true.
  else
     usehist = .true. ! By default we want to use all histograms
  end if

  ! Check for restraints and set usehist to false for restraint
  ! histograms.

  ierror = 0
  histcount2: do while (ierror==0)
    read(iexp, '(a)', iostat=ierror) buffer

    ! Ignore any histograms not set to refine in GSAS
    if (index(buffer, 'HTYP') > 0) then
       buffer = adjustl(buffer(index(buffer, 'HTYP')+4:))
       read(buffer(1:1), '(i1)') ihtyp
       buffer = adjustl(buffer(2:))
       do i = (ihtyp - 1)*12 + 1, min(nhist, ihtyp*12)
          n = index(buffer, ' ')
          if (index(buffer(1:n), '*') > 0) usehist(i) = .false.
          hist_type(i) = buffer(1:3)
          buffer = adjustl(buffer(n+1:))
       end do
    end if

    ! Ignore any restraint histograms
    if (index(buffer, 'HNAM') > 0) then
      if (index(buffer, 'restraints') > 0) then
         read(buffer, *) shortbuffer, ihist
         usehist(ihist) = .false.
      end if
    end if
  end do histcount2
!  rewind(iexp)

  if (ldiag) then
     write(main_output,'(i0,a)') nhist,' histograms of data identified within the EXP file'
     write(main_output,'(a)',advance='no') 'Histograms marked to be processed: '
     do i=1,nhist
        if (usehist(i)) write(main_output, '(i3)', advance='no') i
     end do
     write(main_output,'(a)') ''
  end if

  allocate(cprofile(nhist,0:4))
  allocate(cback(nhist,0:1))
  allocate(nback(nhist))
  allocate(cicons(nhist))
  allocate(chscale(nhist))
  allocate(cbnkpar(nhist))
  allocate(ctrnge(nhist))
  allocate(bank(nhist))
  allocate(nbraggpts(nhist))
  allocate(cabspar(nhist))

! Extract instrument parameters per histogram
  instloop: do ihist = 1,nhist
     if (.not. usehist(ihist)) cycle
!    write(chap,'(i0)') ihist
!    chap = adjustl(chap)
!    chap = 'HAP1 '//trim(chap)//'PRCF'
    write(chap,'("HAP1", i2, "PRCF")') ihist
    ierror = 0
    rewind(iexp)
    readloop1: do while (ierror==0)
    read(iexp,'(a)',iostat=ierror) buffer
    if (ierror/=0) exit readloop1
    if (index(buffer,trim(chap))>0) then
      cprofile(ihist,0) = trim(buffer(13:))
      do i = 1,4
        read(iexp,'(a)') buffer
        cprofile(ihist,i) = trim(buffer(13:))
      end do
!      rewind(iexp)
      exit readloop1
    end if
    end do readloop1
  end do instloop

!  do ihist = 1,nhist
!    do i = 0,4
!    write(6,'(a)') trim(cprofile(ihist,i))
!    end do
!  end do

! Extract bank number for each histogram
  bankloop1: do ihist = 1,nhist
     if (.not. usehist(ihist)) cycle
!    write(cbk,'(i0)') ihist
!    cbk = 'HST  '//trim(cbk)//' BANK'
    write(cbk,'("HST", i3, " BANK")') ihist
    ierror = 0
    rewind(iexp)
    readloop7: do while (ierror==0)
    read(iexp,'(a)',iostat=ierror) buffer
    if (ierror/=0) exit readloop7
    if (index(buffer,trim(cbk))>0) then
      read(buffer(13:),*) bank(ihist)
!      rewind(iexp)
      exit readloop7
    end if
    end do readloop7
  end do bankloop1
  if (ldiag) then
    do ihist = 1,nhist
     if (.not. usehist(ihist)) cycle
      write(main_output,'(a,i0,a,i0)') '... Histogram number ',ihist,' corresponds to data bank number ',bank(ihist)
    end do
  end if

! Extract number of points for each histogram
  bankloop2: do ihist = 1,nhist
     if (.not. usehist(ihist)) cycle
!    write(cpts,'(i0)') ihist
!    cpts = 'HST  '//trim(cpts)//' CHANS'
    write(cpts,'("HST", i3, " CHANS")') ihist
    ierror = 0
    rewind(iexp)
    readloop8: do while (ierror==0)
    read(iexp,'(a)',iostat=ierror) buffer
    if (ierror/=0) exit readloop8
    if (index(buffer,trim(cpts))>0) then
      read(buffer(33:),*) nbraggpts(ihist)
!      rewind(iexp)
      exit readloop8
    end if
    end do readloop8
  end do bankloop2
  if (ldiag) then
    do ihist = 1,nhist
     if (.not. usehist(ihist)) cycle
      write(main_output,'(a,i0,a,i0,a)') '... Histogram number ',ihist,' contains ',nbraggpts(ihist),' points'
    end do
  end if

! Extract background parameters per histogram
  backloop: do ihist = 1,nhist
     if (.not. usehist(ihist)) cycle
!    write(cbk,'(i0)') ihist
!    cbk = 'HST  '//trim(cbk)//'BAKGD'
    write(cbk,'("HST", i3, "BAKGD")') ihist
    ierror = 0
    rewind(iexp)
    readloop2: do while (ierror==0)
    read(iexp,'(a)',iostat=ierror) buffer
    if (ierror/=0) exit readloop2
    if (index(buffer,trim(cbk))>0) then
      cback(ihist,0) = trim(buffer(13:))
      read(cback(ihist,0),*) junk,nback(ihist)
      n = nback(ihist)/4
      if (4*n<nback(ihist)) n = n + 1
      cback(ihist,1) = ''
      do j = 1,n
        read(iexp,'(a)') buffer
        cback(ihist,1) = trim(cback(ihist,1))//' '//trim(buffer(13:))
      end do
!      rewind(iexp)
      exit readloop2
    end if
    end do readloop2
  end do backloop

!  do ihist = 1,nhist
!    write(6,'(i0)') nback(ihist)
!    write(6,'(a)') trim(cback(ihist,1))
!  end do

! Extact ICONS value
  iconsloop: do ihist = 1,nhist
     if (.not. usehist(ihist)) cycle
    write(cbk,'("HST", i3, " ICONS")') ihist
!    cbk = 'HST  '//trim(cbk)//' ICONS'
    ierror = 0
    rewind(iexp)
    readloop3: do while (ierror==0)
    read(iexp,'(a)',iostat=ierror) buffer
    if (ierror/=0) exit readloop3
    if (index(buffer,trim(cbk))>0) then
      cicons(ihist) = trim(buffer(13:))
!          rewind(iexp)
      exit readloop3
    end if
    end do readloop3
  end do iconsloop

!  do ihist = 1,nhist
!    write(6,'(a)') trim(cicons(ihist))
!  end do

! Extact HSCALE value
  hscaleloop: do ihist = 1,nhist
     if (.not. usehist(ihist)) cycle
!    write(cbk,'(i0)') ihist
!    cbk = 'HST  '//trim(cbk)//'HSCALE'
    write(cbk,'("HST", i3, "HSCALE")') ihist
    ierror = 0
    rewind(iexp)
    readloop4: do while (ierror==0)
    read(iexp,'(a)',iostat=ierror) buffer
    if (ierror/=0) exit readloop4
    if (index(buffer,trim(cbk))>0) then
      chscale(ihist) = trim(buffer(13:))
!      rewind(iexp)
      exit readloop4
    end if
    end do readloop4
  end do hscaleloop

!  do ihist = 1,nhist
!    write(6,'(a)') trim(chscale(ihist))
!  end do

! Extact BNKPAR value
  bnkparloop: do ihist = 1,nhist
     if (.not. usehist(ihist)) cycle
!    write(cbk,'(i0)') ihist
!    cbk = 'HST  '//trim(cbk)//'BNKPAR'
    write(cbk,'("HST", i3, "BNKPAR")') ihist
    ierror = 0
    rewind(iexp)
    readloop5: do while (ierror==0)
    read(iexp,'(a)',iostat=ierror) buffer
    if (ierror/=0) exit readloop5
    if (index(buffer,trim(cbk))>0) then
      cbnkpar(ihist) = trim(buffer(13:))
!      rewind(iexp)
      exit readloop5
    end if
    end do readloop5
  end do bnkparloop

  ! Extact ABSCOR value
  absparloop: do ihist = 1,nhist
     if (.not. usehist(ihist)) cycle
    write(cab,'("HST", i3, "ABSCOR")') ihist
    ierror = 0
    rewind(iexp)
    readloop9: do while (ierror==0)
    read(iexp,'(a)',iostat=ierror) buffer
    if (ierror/=0) exit readloop9
    if (index(buffer,trim(cab))>0) then
      cabspar(ihist) = trim(buffer(13:))
      exit readloop9
    end if
    end do readloop9
  end do absparloop

!  do ihist = 1,nhist
!    write(6,'(a)') trim(cbnkpar(ihist))
!  end do

ctrnge = "NULL"

! Extact TRNGE value
  trngeloop: do ihist = 1,nhist
     if (.not. usehist(ihist)) cycle
!    write(cbk,'(i0)') ihist
!    cbk = 'HST  '//trim(cbk)//' TRNGE'
    write(cbk,'("HST", i3, " TRNGE")') ihist
    ierror = 0
    rewind(iexp)
    readloop6: do while (ierror==0)
    read(iexp,'(a)',iostat=ierror) buffer
    if (ierror/=0) exit readloop6
    if (index(buffer,trim(cbk))>0) then
      ctrnge(ihist) = trim(buffer(13:))
!      rewind(iexp)
      exit readloop6
    end if
    end do readloop6
    if (ctrnge(ihist) == "NULL") then
       ! can't use this histogram!
       usehist(ihist) = .false.
       if (ldiag) write (main_output, '("Can''t find TRNGE for histogram ", &
           & i3, ": perhaps it wasn''t used in refinement?")') ihist
    end if
  end do trngeloop

!  do ihist = 1,nhist
!    write(6,'(a)') trim(ctrnge(ihist))
!  end do

! Write the .inst file
  open(iinst,file=trim(fileinst),form='formatted',recl=200,status='unknown')
  write(iinst,'(i0)') count(usehist)
  do ihist = 1,nhist
     if (.not. usehist(ihist)) cycle
    read(cprofile(ihist,0),*) jtype
    read(cicons(ihist),'(3f10.0)') difc,difa,zero
    read(cbnkpar(ihist),*) zzz,theta
    do i = 1,3
    read(cprofile(ihist,i),*) (peak(j+4*(i-1)),j=1,4)
    end do
    write(iinst,'(i0)') ihist
!WS     
    if (hist_type(ihist).eq.'PNT') then 
        read(cbnkpar(ihist),*) zzz,theta
        buff_temp = cicons(ihist)
        read(buff_temp(1:10),'(f10.0)') difc
        read(buff_temp(11:20),'(f10.0)') difa
        read(buff_temp(21:30),'(f10.0)') zero
    write(iinst,*) difc,difa,zero,theta
    elseif (hist_type(ihist).eq.'PXC') then
        buff_temp = cicons(ihist)
        read(buff_temp(1:10),'(f10.0)') lambda
        read(buff_temp(21:30),'(f10.0)') zero
        read(buff_temp(41:50),'(f10.0)') polar
        read(buff_temp(51:55),'(i5)') ipol
        write(iinst,'(3(e12.6e2,2x),i5)') lambda, zero, polar, ipol  
    endif
    if (jtype==1) then  !! GSAS profile type 1
    write(iinst,'(4(e12.6e2,2x))') peak(1),peak(2),peak(3),peak(4)
    write(iinst,'(3(e12.6e2,2x))') peak(5),peak(6),peak(7)
    write(iinst,'(1(e12.6e2,2x))') peak(8)
    else if (jtype==2) then  !! GSAS profile type 2
    write(iinst,'(4(e12.6e2,2x))') peak(1),peak(2),peak(3),peak(4)
    write(iinst,'(3(e12.6e2,2x))') peak(5),peak(6),peak(7)
    write(iinst,'(5(e12.6e2,2x))') peak(8),peak(9),peak(10),peak(11),peak(12)
    else if (jtype==3) then  !! GSAS profile type 3
    write(iinst,'(3(e12.6e2,2x))') peak(1),peak(2),peak(3)
    write(iinst,'(3(e12.6e2,2x))') peak(4),peak(5),peak(6)
    write(iinst,'(5(e12.6e2,2x))') peak(7),peak(8),peak(9),peak(10),peak(11)
    end if
!WS add absorption correction
    buff_temp = cabspar(ihist)
        read(buff_temp(1:15),*) abs_corr
        read(buff_temp(40:45),*) abs_type

    if ((abs_corr.ne.0.0).and.(abs_type.eq.0)) then
        write(iinst,'(1(e12.6e2,2x))') abs_corr
    elseif ((abs_corr.ne.0.0).and.(abs_type.ne.0)) then
      write(*,*) 'You are using GSAS absorption correction type ', abs_type
      write(*,*) 'Sorry this gsas absorption correction type is not supported yet' 
      stop
    endif
  end do
  close(iinst)

!! Ask for the data bank to be used
!  shortbuffer = ''
!  do while ((shortbuffer(1:1)/='H').and.(shortbuffer(1:1)/='h').and. &
!            (shortbuffer(1:1)/='B').and.(shortbuffer(1:1)/='b'))
!    write(6,'(a)', advance="no") 'Do you want to work with Histogram (h/H) or Bank (b/B) numbers: '
!    read(5,'(a)') shortbuffer
!    shortbuffer = adjustl(shortbuffer)
!  end do
!  if ((shortbuffer(1:1)=='H').or.(shortbuffer(1:1)=='h')) then
!    write(6,'(a)', advance="no") 'Please give Histogram number for Bragg profile: '
!    read(5,*) nhistogram
!    if (nhistogram>nhist) then
!      write(6,'(a)') 'You have given too large a histogram number'
!      stop
!    end if
!    nbank = bank(nhistogram)
!  else
!    write(6,'(a)', advance="no") 'Please give Bank number for Bragg profile: '
!    read(5,*) nbank
!    nhistogram = 0
!    do ihist = 1,nhist
!      if (bank(ihist)==nbank) then
!        nhistogram = ihist
!        exit
!      end if
!    end do
!    if (nhistogram==0) then
!      write(6,'(a)') 'You have given a bank number that does not exist'
!      stop
!    end if
!  end if

  n = index(fileroot,'.')
  filelst = filename(1:n)//'LST'
  inquire(file=trim(filelst),exist=l1) ! error fixed AEP, v1.20
  filelst = filename(1:n)//'lst'
  inquire(file=trim(filelst),exist=l2) ! error fixed AEP, v1.20
  if (l1) filelst = trim(filename(1:n)//'LST')
  if (l2) filelst = trim(filename(1:n)//'lst')
  lst_exists = (l1.or.l2)
  if (l1.or.l2) open(ilst,file=trim(filelst),form='formatted',status='old')


! Write the .back files
  do ihist = 1,nhist
     if (.not. usehist(ihist)) cycle
    n = index(fileroot,'.')
    write(fileback,'(i0)') ihist ; fileback = adjustl(fileback)
    fileback = trim(wfolder)//trim(fileroot(1:n-1))//'_'//trim(fileback)//'.back'
    open(iback,file=trim(fileback),form='formatted',status='unknown')
    allocate(back(nback(ihist)))
    read(cback(ihist,1),*) back
    write(iback,'(i0)') nback(ihist)
    do i = 1,nback(ihist)
      write(iback,*) back(i)
    end do
    close(iback)
    deallocate(back)
  end do

  close(iexp)

  if (.not.lst_exists) return

!! Now create the .bragg file from the histogram data contained within the .LST file
!  do ihist = 1,nhist
!    n = index(fileroot,'.')
!    write(filebragg,'(i0)') ihist ; filebragg = adjustl(filebragg)
!    filebragg = trim(wfolder)//trim(fileroot(1:n-1))//'_'//trim(filebragg)//'.bragg'
!    open(ibragg,file=trim(filebragg),form='formatted',status='unknown')
!    rewind(ilst)
!    ierror = 0
!    lstloop1: do while (ierror==0)  ! Position the file just before the first HSTDMP word
!      read(ilst,'(a)',iostat=ierror) buffer
!      if (ierror/=0) exit lstloop1
!      if (index(buffer,'HSTDMP')/=0) then
!        backspace(ilst)
!        exit lstloop1
!      end if
!    end do lstloop1
!  end do

! Check there is actually a histogram in the .LST file

  ierror = 0
  ncheck = 0
  check_loop: do while (ierror==0)
    read(ilst,'(a)',iostat=ierror) buffer
    if (ierror/=0) exit check_loop
    if (index(buffer,'RecNo Code')>0) ncheck = ncheck + 1
  end do check_loop
  rewind(ilst)
  if (ncheck==0) then
    write(6,'(a)') 'No histograms in the LST file'
    return
  end if

  ierror = 0
  lstloop1: do while (ierror==0)  ! Position the file just before the first HSTDMP word
    read(ilst,'(a)',iostat=ierror) buffer
    if (ierror/=0) exit lstloop1
    if (index(buffer,'HSTDMP')/=0) then
      backspace(ilst)
      exit lstloop1
    end if
  end do lstloop1

! Count the total number of points
  npts = 0
  ierror = 0
  lstloop2: do while (ierror==0) 
    read(ilst,'(a)',iostat=ierror) buffer
    if (ierror/=0) exit lstloop2
    if (len_trim(buffer)==0) cycle lstloop2
    if (index(buffer,'HSTDMP')/=0) then
      do while (index(buffer,'RecNo')==0)
        read(ilst,'(a)',iostat=ierror) buffer
        if (ierror/=0) exit lstloop2
      end do
      cycle lstloop2
    end if
    npts = npts + 1
!    read(buffer,*,iostat=ierror2) j
!    if (ierror/=0) exit lstloop2
!    if (j/=npts) then
!      write(6,'(a,i0,1x,i0)') 'Error reading data points ',j,npts
!    end if
  end do lstloop2
  if (ldiag) write(main_output,'(i0,a)') npts,' data points identified from LST file'
  rewind(ilst)

  allocate(tof(npts),ydata(npts),ipt(npts))
  ierror = 0
  lstloop3: do while (ierror==0)
    read(ilst,'(a)',iostat=ierror) buffer
    if (ierror/=0) exit lstloop3
    if (index(buffer,'HSTDMP')/=0) then
      backspace(ilst)
      exit lstloop3
    end if
  end do lstloop3

  j = 0
  ierror = 0
  lstloop4: do while (ierror==0)  ! Position the file just before the first HSTDMP word
    read(ilst,'(a)',iostat=ierror) buffer
    if (ierror/=0) exit lstloop4
    if (len_trim(buffer)==0) cycle lstloop4
    if (index(buffer,'HSTDMP')/=0) then
      do while (index(buffer,'RecNo')==0)
        read(ilst,'(a)',iostat=ierror) buffer
        if (ierror/=0) exit lstloop4
      end do
      cycle lstloop4
    end if
    j = j + 1
    read(buffer,*,iostat=ierror) ipt(j)
    if (ierror/=0) exit lstloop4
    read(buffer(10:),*,iostat=ierror) tof(j),ydata(j)
    if (ierror/=0) exit lstloop4
  end do lstloop4

! Identify the histograms based on the number of points
  cpts = '' ; nhist2 = 0
  do i = 1,npts-1
    if (ipt(i+1)<ipt(i)) then
      n = len_trim(cpts) + 2
      write(cpts(n:),'(i0)') ipt(i)
      nhist2 = nhist2 + 1
    end if
  end do
  nhist2 = nhist2 + 1
  n = len_trim(cpts) + 2
  write(cpts(n:),'(i0)') ipt(npts)

! Need to beware, I am assuming no accidental duplication of histograms in the .lst file
! I wonder whether I can count the number of histograms in the file
  allocate(histpts(nhist2),bank2(nhist2),hist2(nhist2))
  hist2 = -1 ! code that we haven't found this histogram
  read(cpts,*) histpts
  do i = 1,nhist2
    do j = 1,nhist
      if (histpts(i)==nbraggpts(j)) then
        bank2(i) = bank(j)
        hist2(i) = j
        exit
        end if
      end do
    end do

! Here ihist is the actual histogram number, and ihist2 reflects the histogram order
! found in the .LST file
  n2 = 0
! write(6,*) nhist2
  do ihist2=1,nhist2
     if (hist2(ihist2) < 0) cycle
     if (.not. usehist(hist2(ihist2))) cycle
    n1 = n2 + 1
    n2 = n1 + histpts(ihist2) - 1
  n = index(fileroot,'.')
    write(filebragg,'(i0)') hist2(ihist2) ; filebragg = adjustl(filebragg)
  filebragg = trim(wfolder)//trim(fileroot(1:n-1))//'_'//trim(filebragg)//'.bragg'
  open(ibragg,file=trim(filebragg),form='formatted',status='unknown')
    allocate(tofh(histpts(ihist2)),ydatah(histpts(ihist2)))
    tofh = tof(n1:n2) ; ydatah = ydata(n1:n2)
    read(ctrnge(hist2(ihist2)),*) tmin,tmax
    nptsused = histpts(ihist2)
    do i = 1,histpts(ihist2)
      if (tofh(i)<tmin) nptsused = nptsused - 1
      if (tofh(i)>tmax) nptsused = nptsused - 1
    end do
    if (ldiag) write(main_output,'(i0,a,i0)') nptsused,' data points used for bank ',bank2(ihist2)
    read(chscale(hist2(ihist2)),*) scale
  volume = a*b*c*sin(pi*alpha/180.0d0)*sin(pi*beta/180.0d0)*sin(pi*gamma/180.0d0)
    write(ibragg,*) nptsused,bank2(ihist2),scale,volume
!    write(6,'(a)', advance="no") 'Please gave name of data set: '
!    read(5,'(a)') dataname
  dataname = trim(adjustl(metadata_material))
    write(shortbuffer,'(i0)') bank2(ihist2) ; shortbuffer = adjustl(shortbuffer)
  if (len_trim(dataname)>0) then
    dataname = trim(dataname)//', Histogram = '//trim(shortbuffer)
  else
    dataname = 'Histogram = '//trim(shortbuffer)
  end if
    write(shortbuffer,'(i0)') bank2(ihist2) ; shortbuffer = adjustl(shortbuffer)
  dataname = trim(dataname)//', Bank = '//trim(shortbuffer)
  write(ibragg,'(a)') 'Time,   I(obs),  '//trim(dataname)
    do i = 1,histpts(ihist2)
      if (tofh(i)<tmin) cycle
      if (tofh(i)>tmax) cycle
      write(ibragg,*) tofh(i),ydatah(i)
    end do
    close(ibragg)
    deallocate(tofh,ydatah)
  end do

  deallocate(tof,ydata)
  close(ilst)

  return

  end subroutine get_expfile_parameters

!===============================================================================

  subroutine help_info
! ====================

  use version_data

  implicit none

  character(len=80) :: cline

!  write(6,'(/a)') char(27)//'[1m'//'Data2config, version '//trim(version)//char(27)//'[0m'
  write(6,'(a)') '__________________________________________________________________________________________________'

  write(6,'(/a)') char(27)//'[1m'//'Introduction'//char(27)//'[0m'
  write(6,'(a)') 'Data2config is a data conversion utility that was originally developed to'
  write(6,'(a)') 'enable generation of RMCProfile configuration files, but like all these things'
  write(6,'(a)') 'it has been extended, particularly to generate configuration files from CIF'
  write(6,'(a/)') 'format, and to help create file for DLPOLY.'
  write(6,'(a)') 'Like many utilities, the data2config program uses flags to control its operation,'
  write(6,'(a)') 'some of which have arguments, and has one main argument that is the file to be'
  write(6,'(a)') 'converted into a new configuration format'
  write(6,'(a)') '__________________________________________________________________________________________________'

  write(6,'(/a)', advance="no") '... Please press return to continue, Q or q to quit: '
  read(5,'(a)') cline
  if (len_trim(cline)>0) then
    if (index(cline,'q')>0) stop
    if (index(cline,'Q')>0) stop
    write(6,'(a/)') '... Response not recognised so continuing'
  end if

  write(6,'(/a)') char(27)//'[1m'//'Input data formats'//char(27)//'[0m'
  write(6,'(a)') 'Initial configuration file formats are understood by their tag. Recognised tags'
  write(6,'(a)') 'are:'
  write(6,'(a)') '   .TBL or .tbl   - table files produced by GSAS'
  write(6,'(a)') '   .CIF or .cif   - crystal information format'
  write(6,'(a)') '   .CELL or .cell - Discus structure frmat'
  write(6,'(a)') '   .SFF or .sff   - our own simple file format'
  write(6,'(a)') '   .CFG or .cfg   - RMC v3 configuration file'
  write(6,'(a)') '   .HIS or .his   - RMC v3 histogram file'
  write(6,'(a)') '   .CSSR or .cssr - CSSR file'
  write(6,'(a)') '   .RES or .res - SHELX/CASTEP res file'
!  write(6,'(a)') '   .CML or .cml   - Chemical Markup Language file *'
!  write(6,'(a)') '   .XML or .xml   - assumed CML language *'
  write(6,'(a)') '__________________________________________________________________________________________________'
!  write(6,'(a)') '* denotes this functionality is not yet implemented but will be soon'

  write(6,'(/a)', advance="no") '... Please press return to continue, Q or q to quit: '
  read(5,'(a)') cline
  if (len_trim(cline)>0) then
    if (index(cline,'q')>0) stop
    if (index(cline,'Q')>0) stop
    write(6,'(a/)') '... Response not recognised so continuing'
  end if

  write(6,'(/a)') char(27)//'[1m'//'Configuration format flags'//char(27)//'[0m'
  write(6,'(a)') '-atomeye     write atomeye configuration file'
!  write(6,'(a)') '-cml         write CML file *'
  write(6,'(a)') '-crystal     write crystal input file'
  write(6,'(a)') '-cif         write cif output file for visualisation'
  write(6,'(a)') '-cmtx        write CrystalMaker text output file for visualisation *'
  write(6,'(a)') '-cssr        write cssr output file for visualisation'
  write(6,'(a)') '-dlpoly      write DLPOLY CONFIG configuration file'
  write(6,'(a)') '-gasp        write GASP configuration file (v3)'
  write(6,'(a)') '-gulpc       write GULP cartesian coordinates file'
  write(6,'(a)') '-gulpf       write GULP fractional coordinates file'
  write(6,'(a)') '-rmc3        write RMCProfile v3 configuration file (extension .cfg)'
  write(6,'(a)') '-rmc6f       write RMCProfile v6f configuration file'
  write(6,'(a)') '-rmc6o       write RMCProfile v6o configuration file'
  write(6,'(a)') '-www         write file in WWW format for use with PAFmaker'  
  write(6,'(a)') '-xtl         write XTL-format configuration file'
  write(6,'(a)') 'Note that you are not able to select both -atomeye and -rmc6f flags at the'
  write(6,'(a)') 'same time, because both use the same filename extension. The simple solution'
  write(6,'(a)') 'is to run data2config twice.'
  write(6,'(a)') '* Note that if a magnetic .cfg file is present, the -cmtx option will allow'
  write(6,'(a)') '  for generation of a set of vectors corresponding to the spins'
  write(6,'(a)') '__________________________________________________________________________________________________'

  write(6,'(/a)', advance="no") '... Please press return to continue, Q or q to quit: '
  read(5,'(a)') cline
  if (len_trim(cline)>0) then
    if (index(cline,'q')>0) stop
    if (index(cline,'Q')>0) stop
    write(6,'(a/)') '... Response not recognised so continuing'
  end if

  write(6,'(/a)') char(27)//'[1m'//'Structure manipulation flags'//char(27)//'[0m'
  
  write(6,'(a)') '-alighn      Align the structure compared to a reference file, requires a'
  write(6,'(a)') '             reference rmc6f file input'  
  write(6,'(a)') '-analysis    Generate simple bond analysis, requires file input'
!  write(6,'(a)') '-cml         write CML file *'
  write(6,'(a)') '-diag        write diagnostics file'
  write(6,'(a)') '-gridcell    give a new size of the grid cell used internally to speed up'
  write(6,'(a)') '             calculations of neighbour lists. The default value is 5.0 Å'
  write(6,'(a)') '             (requires argument for the grid cell size)'
  write(6,'(a)') '-nano        produce a spherical nanoparticle of a required radius'
  write(6,'(a)') '             (requires argument for the radius)'
  write(6,'(a)') '-nanobox     Specify the multipliers for the nanoparticle box.'
  write(6,'(a)') '             If not provided, [4 4 4] will be used.'
  write(6,'(a)') '-nanoscale   scale a nanoparticle size and move the ligands radially outwards'
  write(6,'(a)') '-noannotate  request no annotations to files'
  write(6,'(a)') '-order       order the atoms by user-supplied list within [...] brackets *'
  write(6,'(a)') '-origin      reset the orgin so that the centre is really the centre'
  write(6,'(a)') '-orient      reorient molecules according to a supplied argument:'
  write(6,'(a)') '             random = by a random rotation'
  write(6,'(a)') '             orth or 90 = by a random rotation of 90 degrees'
  write(6,'(a)') '             flipx or 180x = by a random rotation of 180 degrees about the x axis'
  write(6,'(a)') '             flipy or 180y = by a random rotation of 180 degrees about the y axis'
  write(6,'(a)') '             flipz or 180z = by a random rotation of 180 degrees about the z axis'
  write(6,'(a)') '             flipr or 180r = by a random rotation of 180 degrees about a random axis'
  write(6,'(a)') '-rect        convert hexagonal cell to a C-centred orthrhombic cell'
  write(6,'(a)') '-reduce      move all atoms within a supercell configuration back into the'
  write(6,'(a)') '             original unit cell. In this case, the -supercell keyword gives'
  write(6,'(a)') '             the supercell integers for scaling down rather than up'
  write(6,'(a)') '-reorder     Reorder the atoms by user-supplied list within [...] brackets *'
  write(6,'(a)') '             This option enforces -one, that is, no supercell is created.'
  write(6,'(a)') '-replace     option to replace selected atoms with another type of atom:'
  write(6,'(a)') '             give atom to replace, fraction to be replaced, and new atom'
  write(6,'(a)') '             within [...] brackets, eg as [Na 0.5 K]. You can give labels'
  write(6,'(a)') '             with the first atom, eg Na2. You can give any number of replacement'
  write(6,'(a)') '             separated by colons.'
  write(6,'(a)') '-rsame       minimum separation of distinct atoms (default = 0.01 Angstrom)'
  write(6,'(a)') '-size        deduce supercell multipliers from approximate required cell size'
  write(6,'(a)') '             (requires argument for the approximate size)'
  write(6,'(a)') '-shift       Shifts all fractional coordinates to between 0 and +1' 
  write(6,'(a)') '-shifts      Shifts all coordinates by an amount given as a parameter (in Angstron)' 
  write(6,'(a)') '             Also -shifty and -shiftz'  
  write(6,'(a)') '-sort        sort the atoms by atomic number'
  write(6,'(a)') '-supercell   option to provide supercell vector as 3 integers enclosed'
  write(6,'(a)') '             within [....] brackets*'
!  write(6,'(a)') '-transform   option to provide transformation matrix as 9 integers enclosed'
!  write(6,'(a)') '             within [....] brackets row-wise'
  write(6,'(a)') '-Uiso        apply Gaussian distribution displacement for all atoms'
  write(6,'(a)') '             The value of Uiso for each atom type should be given in [...]'
  write(6,'(a)') '-Biso        apply Gaussian distribution displacement for all atoms'
  write(6,'(a)') '             The value of Biso for each atom type should be given in [...]'  
  write(6,'(a)') '-vacancy     option to replace some atoms by vacant sites at random; give atom'
  write(6,'(a)') '             label with optional site number attached, eg Na2, followed by fraction'
  write(6,'(a)') '             of vacancies to create, in [...] brackets; you can give any number of'
  write(6,'(a)') '             atom/fraction pairs you like, separating pairs with colons (note that'
  write(6,'(a)') '             the labels give the sequencenumber of each atom in the original input'
  write(6,'(a)') '             file, eg Na1, Na2 for two Na atoms, followed by O1, O2, O3 for three'
  write(6,'(a)') '             oxygen atoms)'
  write(6,'(a)') '-delete      when used with -vacancy, this deletes all vacant sites'
  write(6,'(a)') '-bank        give the bank number of the GSAS histogram'
  write(6,'(a)') 'If the input file is of type .TBL or .tbl as generated by GSAS, there are two'
  write(6,'(a)') 'additional possibilities:'
  write(6,'(a)') '1. If the GSAS .EXP or .exp file is present, the .back and .inst files will be'
  write(6,'(a)') '   created automatically.'
  write(6,'(a)') '2. If the GSAS .LST or .lst file is present and the data histogram has been'
  write(6,'(a)') '   written to it, the .bragg file will be created automatically.'
  write(6,'(a)') '* Footnote, if parameters are not supplied within the [...] brackets, they can be'
  write(6,'(a)') 'entered interactively in response to questions'
  write(6,'(a)') '__________________________________________________________________________________________________'

  write(6,'(/a)', advance="no") '... Please press return to continue, Q or q to quit: '
  read(5,'(a)') cline
  if (len_trim(cline)>0) then
    if (index(cline,'q')>0) stop
    if (index(cline,'Q')>0) stop
    write(6,'(a/)') '... Response not recognised so continuing'
  end if

  write(6,'(/a)') char(27)//'[1m'//'DLPOLY flags'//char(27)//'[0m'
  write(6,'(a)') '-cssr        write cssr file'
  write(6,'(a)') '-dlpoly      write DLPOLY CONFIG file'
  write(6,'(a)') '-diag        write diagnostics file'
  write(6,'(a)') '-field       write intramolecular part of DLPOLY FIELD file, requires file input [1]'
  write(6,'(a)') '-order       order the atoms by user-supplied list within [...] brackets [2]'
  write(6,'(a)') '-sort        sort the atoms by atomic number'
  write(6,'(a)') '[1] If the -field option is used, this requires a file name as an argument.'
  write(6,'(a)') '    This file contains the following optional lines:'
  write(6,'(a)') '    rigid <2 x atoms> <distance> <variation>'
  write(6,'(a)') '    bond <2 x atoms> <distance> <variation> : <DLPOLY potential>'
  write(6,'(a)') '    angle <3 x atoms> <2 x distances> <variation> <angle> <variation> : <DLPOLY potential>'
  write(6,'(a)') '    torsion <4 x atoms> <4 x distances> <variation> <2 x angles> <variation> : <DLPOLY potential>'
  write(6,'(a)') '    inversion <4 x atoms> <3 x distances> <variation> <2 x angles> <variation> : <DLPOLY potential>'
  write(6,'(a)') '    The potentials part is the statement required in the DLPOLY field file.'
  write(6,'(a)') '[2] If parameters are not supplied within the [...] brackets, they can be'
  write(6,'(a)') '    entered interactively in response to questions'
  write(6,'(a)') '__________________________________________________________________________________________________'

  write(6,'(/a)', advance="no") '... Please press return to continue, Q or q to quit: '
  read(5,'(a)') cline
  if (len_trim(cline)>0) then
    if (index(cline,'q')>0) stop
    if (index(cline,'Q')>0) stop
    write(6,'(a/)') '... Response not recognised so continuing'
  end if

  write(6,'(/a)') char(27)//'[1m'//'Analysis flags'//char(27)//'[0m'
  write(6,'(a)') '-analysis    Generate simple bond analysis, requires file input.'
  write(6,'(a)') '             The file contains one line per atom pair, with data'
  write(6,'(a)') '             <2 x atoms> <min distance> <max distance>'
  write(6,'(a)') '             This computes a histogram of coordination numbers.'
  write(6,'(a)') '-cmbonds     Generate list of interatomic distances from a pair list generated by CrystalMaker.'
  write(6,'(a)') '             Requires CrystalMaker text file as argument'
  write(6,'(a)') '-compare     Compare the input rmc6f file with another rmc6f file, which is provided as'
  write(6,'(a)') '             the supplied argument. Note that this does not check for consistency of cell'
  write(6,'(a)') '             size, but it does check on the number of atoms.'
  write(6,'(a)') '-dlanal      Analysis of DL_POLY CONFIG or REVCON file based on the bonds and angles given.'
  write(6,'(a)') '             in the FIELD file.'
  write(6,'(a)') '-listf       Give a file containing a list of atoms required for the origins of the analysis'
  write(6,'(a)') '-limits      Define the range of x, y and z of the target atoms in calculatins of the PDF or'
  write(6,'(a)') '             the bond analysis. The limits are given in [...] brackets in terms of fractional'
  write(6,'(a)') '             coordinates, in the order of xmin, xmax, ymin, ymax, zmin and zmax, where min and'
  write(6,'(a)') '             max define the minimum and maximum values of x, y or z'
  write(6,'(a)') '-molecules   locates molecules in a configuration, requires file input [1]'
  write(6,'(a)') '-onecell     Use only the atoms in one unit cell for analysis'
  write(6,'(a)') '-pdf         Calculates the pair distribution function and corresponding scattering functions'
  write(6,'(a)') '             for both neutron and x-ray scattering. The maximum distance is supplied as argument.'
  write(6,'(a)') '             This facility works in three cases; 1) for a crystal entry using -one as a keyword;'
  write(6,'(a)') '             2) when the -nano flag is used; and 3) when the input structure is recognised as a'
  write(6,'(a)') '             configuration (RMC or MD).'
  write(6,'(a)') '-pdfwidth    Give a value for the width of the PDF if automatic broadening is required, supplying'
  write(6,'(a)') '             as the following argument.'
  write(6,'(a)') '-rmcanal     Performs analysis of the bond and bond angle distributions as defined by the molecule'
  write(6,'(a)') '             .bonds and .triplets files'
  write(6,'(a)') '-ylm         Calculates Kubic harmonic coefficients and orientational variables for molecules'
  write(6,'(a)') '             on sites of cubic symmetry. This requires use of the -molecules flag if the file'
  write(6,'(a)') '             the file molecules.dat does not exist. It generates a file ylm.out' 
  write(6,'(a)') '[1] This file contains the following optional lines:'
  write(6,'(a)') '   bond <2 x atoms> <distance> <variation>'
  write(6,'(a)') '__________________________________________________________________________________________________'

  write(6,'(/)')

  stop
  end subroutine help_info

! >>>>>>>>>>>>>>>>>>>> Yuanpeng here >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! This subroutine determines whether the given string is numeric. 
! The code was copied (available under the GNU free license) from the following website:
! https://rosettacode.org/wiki/Determine_if_a_string_is_numeric#Fortran
   FUNCTION is_numeric(string)
     IMPLICIT NONE
     CHARACTER(len=*), INTENT(IN) :: string
     LOGICAL :: is_numeric
     REAL :: x
     INTEGER :: e
     READ(string,*,IOSTAT=e) x
     is_numeric = e == 0
   END FUNCTION is_numeric
! <<<<<<<<<<<<<<<<<<<< Yuanpeng finishes here <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
