# $Id: spatt.spec.in 438 2005-07-20 12:09:50Z mhoebeke $
%define name spatt
%define version 1.2.2
%define release 1mdk

Summary: A set of programs for performing statistics on words or  patterns.
Name: %{name}
Version: %{version}
Release: %{release}
Source0: http://stat.genopole.cnrs.fr/%{name}/%{name}-%{version}.tar.bz2
License: GPL
Group: Sciences/Biology
Url: http://stat.genopole.cnrs.fr/%{name}
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-buildroot
BuildRequires: gsl-devel gcc gcc-c++ gcc-g77

%description
SPATT provides a set of programs performing statistical analysis of
words or patterns in biological sequences.


%prep

%setup -q

%build

%configure

%make

%install
rm -rf $RPM_BUILD_ROOT
%makeinstall

%clean
rm -rf $RPM_BUILD_ROOT

%files 
%defattr(-,root,root)
%{_bindir}/cpspatt
%{_bindir}/gspatt
%{_bindir}/ldspatt
%{_bindir}/sspatt
%{_bindir}/xspatt
%{_mandir}/man1/cpspatt.1.bz2
%{_mandir}/man1/gspatt.1.bz2
%{_mandir}/man1/ldspatt.1.bz2
%{_mandir}/man1/spatt.1.bz2
%{_mandir}/man1/sspatt.1.bz2
%{_mandir}/man1/xspatt.1.bz2
%doc AUTHORS COPYING README NEWS

%changelog
* Tue Jan 18 2005 Hoebeke Mark <Mark.Hoebeke@jouy.inra.fr> 1.2.2-1mdk
- First attempt to build SPATT package.


