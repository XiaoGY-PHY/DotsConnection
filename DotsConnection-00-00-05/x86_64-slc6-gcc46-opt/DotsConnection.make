#-- start of make_header -----------------

#====================================
#  Library DotsConnection
#
#   Generated Thu Jul  8 16:09:03 2021  by mg20220135
#
#====================================

include ${CMTROOT}/src/Makefile.core

ifdef tag
CMTEXTRATAGS = $(tag)
else
tag       = $(CMTCONFIG)
endif

cmt_DotsConnection_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_DotsConnection_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_DotsConnection

DotsConnection_tag = $(tag)

#cmt_local_tagfile_DotsConnection = $(DotsConnection_tag)_DotsConnection.make
cmt_local_tagfile_DotsConnection = $(bin)$(DotsConnection_tag)_DotsConnection.make

else

tags      = $(tag),$(CMTEXTRATAGS)

DotsConnection_tag = $(tag)

#cmt_local_tagfile_DotsConnection = $(DotsConnection_tag).make
cmt_local_tagfile_DotsConnection = $(bin)$(DotsConnection_tag).make

endif

include $(cmt_local_tagfile_DotsConnection)
#-include $(cmt_local_tagfile_DotsConnection)

ifdef cmt_DotsConnection_has_target_tag

cmt_final_setup_DotsConnection = $(bin)setup_DotsConnection.make
cmt_dependencies_in_DotsConnection = $(bin)dependencies_DotsConnection.in
#cmt_final_setup_DotsConnection = $(bin)DotsConnection_DotsConnectionsetup.make
cmt_local_DotsConnection_makefile = $(bin)DotsConnection.make

else

cmt_final_setup_DotsConnection = $(bin)setup.make
cmt_dependencies_in_DotsConnection = $(bin)dependencies.in
#cmt_final_setup_DotsConnection = $(bin)DotsConnectionsetup.make
cmt_local_DotsConnection_makefile = $(bin)DotsConnection.make

endif

#cmt_final_setup = $(bin)setup.make
#cmt_final_setup = $(bin)DotsConnectionsetup.make

#DotsConnection :: ;

dirs ::
	@if test ! -r requirements ; then echo "No requirements file" ; fi; \
	  if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi

javadirs ::
	@if test ! -d $(javabin) ; then $(mkdir) -p $(javabin) ; fi

srcdirs ::
	@if test ! -d $(src) ; then $(mkdir) -p $(src) ; fi

help ::
	$(echo) 'DotsConnection'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = DotsConnection/
#DotsConnection::
#	@if test ! -d $(bin)$(binobj) ; then $(mkdir) -p $(bin)$(binobj) ; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)$(binobj)
endif

${CMTROOT}/src/Makefile.core : ;
ifdef use_requirements
$(use_requirements) : ;
endif

#-- end of make_header ------------------
#-- start of libary_header ---------------

DotsConnectionlibname   = $(bin)$(library_prefix)DotsConnection$(library_suffix)
DotsConnectionlib       = $(DotsConnectionlibname).a
DotsConnectionstamp     = $(bin)DotsConnection.stamp
DotsConnectionshstamp   = $(bin)DotsConnection.shstamp

DotsConnection :: dirs  DotsConnectionLIB
	$(echo) "DotsConnection ok"

#-- end of libary_header ----------------

DotsConnectionLIB :: $(DotsConnectionlib) $(DotsConnectionshstamp)
	@/bin/echo "------> DotsConnection : library ok"

$(DotsConnectionlib) :: $(bin)DotsConnection_load.o $(bin)DotsConnection_entries.o $(bin)DotsConnection.o $(bin)TrkFitFun.o $(bin)DotsHelixFitter.o
	$(lib_echo) library
	$(lib_silent) cd $(bin); \
	  $(ar) $(DotsConnectionlib) $?
	$(lib_silent) $(ranlib) $(DotsConnectionlib)
	$(lib_silent) cat /dev/null >$(DotsConnectionstamp)

#------------------------------------------------------------------
#  Future improvement? to empty the object files after
#  storing in the library
#
##	  for f in $?; do \
##	    rm $${f}; touch $${f}; \
##	  done
#------------------------------------------------------------------

$(DotsConnectionlibname).$(shlibsuffix) :: $(DotsConnectionlib) $(DotsConnectionstamps)
	$(lib_silent) cd $(bin); QUIET=$(QUIET); $(make_shlib) "$(tags)" DotsConnection $(DotsConnection_shlibflags)

$(DotsConnectionshstamp) :: $(DotsConnectionlibname).$(shlibsuffix)
	@if test -f $(DotsConnectionlibname).$(shlibsuffix) ; then cat /dev/null >$(DotsConnectionshstamp) ; fi

DotsConnectionclean ::
	$(cleanup_echo) objects
	$(cleanup_silent) cd $(bin); /bin/rm -f $(bin)DotsConnection_load.o $(bin)DotsConnection_entries.o $(bin)DotsConnection.o $(bin)TrkFitFun.o $(bin)DotsHelixFitter.o

#-----------------------------------------------------------------
#
#  New section for automatic installation
#
#-----------------------------------------------------------------

ifeq ($(INSTALLAREA),)
installarea = $(CMTINSTALLAREA)
else
ifeq ($(findstring `,$(INSTALLAREA)),`)
installarea = $(shell $(subst `,, $(INSTALLAREA)))
else
installarea = $(INSTALLAREA)
endif
endif

install_dir = ${installarea}/${CMTCONFIG}/lib
DotsConnectioninstallname = $(library_prefix)DotsConnection$(library_suffix).$(shlibsuffix)

DotsConnection :: DotsConnectioninstall

install :: DotsConnectioninstall

DotsConnectioninstall :: $(install_dir)/$(DotsConnectioninstallname)
	@if test ! "${installarea}" = ""; then\
	  echo "installation done"; \
	fi

$(install_dir)/$(DotsConnectioninstallname) :: $(bin)$(DotsConnectioninstallname)
	@if test ! "${installarea}" = ""; then \
	  cd $(bin); \
	  if test ! "$(install_dir)" = ""; then \
	    if test ! -d "$(install_dir)"; then \
	      mkdir -p $(install_dir); \
	    fi ; \
	    if test -d "$(install_dir)"; then \
	      echo "Installing library $(DotsConnectioninstallname) into $(install_dir)"; \
	      if test -e $(install_dir)/$(DotsConnectioninstallname); then \
	        $(cmt_uninstall_area_command) $(install_dir)/$(DotsConnectioninstallname); \
	        $(cmt_uninstall_area_command) $(install_dir)/$(DotsConnectioninstallname).cmtref; \
	      fi; \
	      $(cmt_install_area_command) `pwd`/$(DotsConnectioninstallname) $(install_dir)/$(DotsConnectioninstallname); \
	      echo `pwd`/$(DotsConnectioninstallname) >$(install_dir)/$(DotsConnectioninstallname).cmtref; \
	    fi \
          else \
	    echo "Cannot install library $(DotsConnectioninstallname), no installation directory specified"; \
	  fi; \
	fi

DotsConnectionclean :: DotsConnectionuninstall

uninstall :: DotsConnectionuninstall

DotsConnectionuninstall ::
	@if test ! "${installarea}" = ""; then \
	  cd $(bin); \
	  if test ! "$(install_dir)" = ""; then \
	    if test -d "$(install_dir)"; then \
	      echo "Removing installed library $(DotsConnectioninstallname) from $(install_dir)"; \
	      $(cmt_uninstall_area_command) $(install_dir)/$(DotsConnectioninstallname); \
	      $(cmt_uninstall_area_command) $(install_dir)/$(DotsConnectioninstallname).cmtref; \
	    fi \
          else \
	    echo "Cannot uninstall library $(DotsConnectioninstallname), no installation directory specified"; \
	  fi; \
	fi




#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),DotsConnectionclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)DotsConnection_load.d

$(bin)$(binobj)DotsConnection_load.d :

$(bin)$(binobj)DotsConnection_load.o : $(cmt_final_setup_DotsConnection)

$(bin)$(binobj)DotsConnection_load.o : $(src)DotsConnection_load.cxx
	$(cpp_echo) $(src)DotsConnection_load.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(DotsConnection_pp_cppflags) $(lib_DotsConnection_pp_cppflags) $(DotsConnection_load_pp_cppflags) $(use_cppflags) $(DotsConnection_cppflags) $(lib_DotsConnection_cppflags) $(DotsConnection_load_cppflags) $(DotsConnection_load_cxx_cppflags)  $(src)DotsConnection_load.cxx
endif
endif

else
$(bin)DotsConnection_dependencies.make : $(DotsConnection_load_cxx_dependencies)

$(bin)DotsConnection_dependencies.make : $(src)DotsConnection_load.cxx

$(bin)$(binobj)DotsConnection_load.o : $(DotsConnection_load_cxx_dependencies)
	$(cpp_echo) $(src)DotsConnection_load.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(DotsConnection_pp_cppflags) $(lib_DotsConnection_pp_cppflags) $(DotsConnection_load_pp_cppflags) $(use_cppflags) $(DotsConnection_cppflags) $(lib_DotsConnection_cppflags) $(DotsConnection_load_cppflags) $(DotsConnection_load_cxx_cppflags)  $(src)DotsConnection_load.cxx

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),DotsConnectionclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)DotsConnection_entries.d

$(bin)$(binobj)DotsConnection_entries.d :

$(bin)$(binobj)DotsConnection_entries.o : $(cmt_final_setup_DotsConnection)

$(bin)$(binobj)DotsConnection_entries.o : $(src)DotsConnection_entries.cxx
	$(cpp_echo) $(src)DotsConnection_entries.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(DotsConnection_pp_cppflags) $(lib_DotsConnection_pp_cppflags) $(DotsConnection_entries_pp_cppflags) $(use_cppflags) $(DotsConnection_cppflags) $(lib_DotsConnection_cppflags) $(DotsConnection_entries_cppflags) $(DotsConnection_entries_cxx_cppflags)  $(src)DotsConnection_entries.cxx
endif
endif

else
$(bin)DotsConnection_dependencies.make : $(DotsConnection_entries_cxx_dependencies)

$(bin)DotsConnection_dependencies.make : $(src)DotsConnection_entries.cxx

$(bin)$(binobj)DotsConnection_entries.o : $(DotsConnection_entries_cxx_dependencies)
	$(cpp_echo) $(src)DotsConnection_entries.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(DotsConnection_pp_cppflags) $(lib_DotsConnection_pp_cppflags) $(DotsConnection_entries_pp_cppflags) $(use_cppflags) $(DotsConnection_cppflags) $(lib_DotsConnection_cppflags) $(DotsConnection_entries_cppflags) $(DotsConnection_entries_cxx_cppflags)  $(src)DotsConnection_entries.cxx

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),DotsConnectionclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)DotsConnection.d

$(bin)$(binobj)DotsConnection.d :

$(bin)$(binobj)DotsConnection.o : $(cmt_final_setup_DotsConnection)

$(bin)$(binobj)DotsConnection.o : $(src)DotsConnection.cxx
	$(cpp_echo) $(src)DotsConnection.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(DotsConnection_pp_cppflags) $(lib_DotsConnection_pp_cppflags) $(DotsConnection_pp_cppflags) $(use_cppflags) $(DotsConnection_cppflags) $(lib_DotsConnection_cppflags) $(DotsConnection_cppflags) $(DotsConnection_cxx_cppflags)  $(src)DotsConnection.cxx
endif
endif

else
$(bin)DotsConnection_dependencies.make : $(DotsConnection_cxx_dependencies)

$(bin)DotsConnection_dependencies.make : $(src)DotsConnection.cxx

$(bin)$(binobj)DotsConnection.o : $(DotsConnection_cxx_dependencies)
	$(cpp_echo) $(src)DotsConnection.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(DotsConnection_pp_cppflags) $(lib_DotsConnection_pp_cppflags) $(DotsConnection_pp_cppflags) $(use_cppflags) $(DotsConnection_cppflags) $(lib_DotsConnection_cppflags) $(DotsConnection_cppflags) $(DotsConnection_cxx_cppflags)  $(src)DotsConnection.cxx

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),DotsConnectionclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)TrkFitFun.d

$(bin)$(binobj)TrkFitFun.d :

$(bin)$(binobj)TrkFitFun.o : $(cmt_final_setup_DotsConnection)

$(bin)$(binobj)TrkFitFun.o : $(src)TrkFitFun.cxx
	$(cpp_echo) $(src)TrkFitFun.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(DotsConnection_pp_cppflags) $(lib_DotsConnection_pp_cppflags) $(TrkFitFun_pp_cppflags) $(use_cppflags) $(DotsConnection_cppflags) $(lib_DotsConnection_cppflags) $(TrkFitFun_cppflags) $(TrkFitFun_cxx_cppflags)  $(src)TrkFitFun.cxx
endif
endif

else
$(bin)DotsConnection_dependencies.make : $(TrkFitFun_cxx_dependencies)

$(bin)DotsConnection_dependencies.make : $(src)TrkFitFun.cxx

$(bin)$(binobj)TrkFitFun.o : $(TrkFitFun_cxx_dependencies)
	$(cpp_echo) $(src)TrkFitFun.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(DotsConnection_pp_cppflags) $(lib_DotsConnection_pp_cppflags) $(TrkFitFun_pp_cppflags) $(use_cppflags) $(DotsConnection_cppflags) $(lib_DotsConnection_cppflags) $(TrkFitFun_cppflags) $(TrkFitFun_cxx_cppflags)  $(src)TrkFitFun.cxx

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),DotsConnectionclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)DotsHelixFitter.d

$(bin)$(binobj)DotsHelixFitter.d :

$(bin)$(binobj)DotsHelixFitter.o : $(cmt_final_setup_DotsConnection)

$(bin)$(binobj)DotsHelixFitter.o : $(src)DotsHelixFitter.cxx
	$(cpp_echo) $(src)DotsHelixFitter.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(DotsConnection_pp_cppflags) $(lib_DotsConnection_pp_cppflags) $(DotsHelixFitter_pp_cppflags) $(use_cppflags) $(DotsConnection_cppflags) $(lib_DotsConnection_cppflags) $(DotsHelixFitter_cppflags) $(DotsHelixFitter_cxx_cppflags)  $(src)DotsHelixFitter.cxx
endif
endif

else
$(bin)DotsConnection_dependencies.make : $(DotsHelixFitter_cxx_dependencies)

$(bin)DotsConnection_dependencies.make : $(src)DotsHelixFitter.cxx

$(bin)$(binobj)DotsHelixFitter.o : $(DotsHelixFitter_cxx_dependencies)
	$(cpp_echo) $(src)DotsHelixFitter.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(DotsConnection_pp_cppflags) $(lib_DotsConnection_pp_cppflags) $(DotsHelixFitter_pp_cppflags) $(use_cppflags) $(DotsConnection_cppflags) $(lib_DotsConnection_cppflags) $(DotsHelixFitter_cppflags) $(DotsHelixFitter_cxx_cppflags)  $(src)DotsHelixFitter.cxx

endif

#-- end of cpp_library ------------------
#-- start of cleanup_header --------------

clean :: DotsConnectionclean ;
#	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(DotsConnection.make) $@: No rule for such target" >&2
else
.DEFAULT::
	$(error PEDANTIC: $@: No rule for such target)
endif

DotsConnectionclean ::
#-- end of cleanup_header ---------------
#-- start of cleanup_library -------------
	$(cleanup_echo) library DotsConnection
	-$(cleanup_silent) cd $(bin); /bin/rm -f $(library_prefix)DotsConnection$(library_suffix).a $(library_prefix)DotsConnection$(library_suffix).s? DotsConnection.stamp DotsConnection.shstamp
#-- end of cleanup_library ---------------
