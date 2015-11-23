Name: fgsl	
Version: 1.0.0
Release:	1%{?dist}
Summary: Fortran 2003 bindings for the GNU Scientific Library	

Group:	Applications/Engineering and Scientific
License: GPL
URL: http://www.lrz.de/services/software/mathematik/gsl/fortran/index.html	
Source0: fgsl-%{version}.tar.gz
BuildRoot:	%(mktemp -ud %{_tmppath}/%{name}-%{version}-%{release}-XXXXXX)

BuildRequires:	gsl-devel >= 1.13 gcc gcc-gfortran glibc-headers glibc-devel
Requires: gsl >= 1.13 libgfortran

%description
Fortran 2003 bindings for the GNU Scientific Library


%package devel
Summary: Development files for FGSL
Requires: gsl-devel >= 1.13 pkgconfig fgsl

%description devel
Development files for FGSL

%package doc
Summary: Documentation for FGSL
Requires: fgsl-devel

%description doc
Documentation and examples for FGSL

%prep
%setup -q


%build
%configure FC=gfortran CC=gcc

sed -i 's|^hardcode_libdir_flag_spec=.*|hardcode_libdir_flag_spec=""|g' libtool
sed -i 's|^runpath_var=LD_RUN_PATH|runpath_var=DIE_RPATH_DIE|g' libtool

make %{?_smp_mflags}




%install
rm -rf $RPM_BUILD_ROOT
make install DESTDIR=$RPM_BUILD_ROOT
libtool --finish $RPM_BUILD_ROOT%{_libdir}
rm -f $RPM_BUILD_ROOT%{_libdir}/*.la



%clean
rm -rf $RPM_BUILD_ROOT

%post -p /sbin/ldconfig

%postun -p /sbin/ldconfig

%files
%defattr(-,root,root,-)
%{_libdir}/libfgsl.so.*

%files devel
%defattr(-,root,root,-)
%{_libdir}/libfgsl.so
%{_libdir}/libfgsl.a
%{_includedir}/fgsl/fgsl.mod
%{_libdir}/pkgconfig/fgsl.pc

%files doc
%defattr(-,root,root,-)
%{_datadir}/doc/fgsl/*
%{_datadir}/examples/fgsl/*




%changelog
* Tue Feb 11 2014 Tom Schoonjans
- added doc package

* Fri Nov 4 2011 Tom Schoonjans
- gsl version requirement

* Mon Apr 18 2011 Tom Schoonjans
- File creation

