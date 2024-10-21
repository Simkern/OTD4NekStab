INC_SRC	  := ./inc_src
FRAME_SRC := ./frame

$(OBJDIR)/frame.o               :$(FRAME_SRC)/frame.f;                  $(F77) -c $(FL2) -I$(INC_SRC) $< -o $@
$(OBJDIR)/mntrlog_block.o       :$(FRAME_SRC)/mntrlog_block.f;          $(F77) -c $(FL2) -I$(INC_SRC) $< -o $@
$(OBJDIR)/mntrlog.o             :$(FRAME_SRC)/mntrlog.f;                $(F77) -c $(FL2) -I$(INC_SRC) $< -o $@
$(OBJDIR)/mntrtmr_block.o       :$(FRAME_SRC)/mntrtmr_block.f;          $(F77) -c $(FL2) -I$(INC_SRC) $< -o $@
$(OBJDIR)/mntrtmr.o             :$(FRAME_SRC)/mntrtmr.f;                $(F77) -c $(FL2) -I$(INC_SRC) $< -o $@
$(OBJDIR)/rprm_block.o          :$(FRAME_SRC)/rprm_block.f;             $(F77) -c $(FL2) -I$(INC_SRC) $< -o $@
$(OBJDIR)/rprm.o                :$(FRAME_SRC)/rprm.f;                   $(F77) -c $(FL2) -I$(INC_SRC) $< -o $@
$(OBJDIR)/io_tools.o            :$(FRAME_SRC)/io_tools.f;               $(F77) -c $(FL2) -I$(INC_SRC) $< -o $@
$(OBJDIR)/io_tools_block.o      :$(FRAME_SRC)/io_tools_block.f;         $(F77) -c $(FL2) -I$(INC_SRC) $< -o $@
$(OBJDIR)/math_tools.o          :math_tools.f;                          $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/otd.o                 :otd.f;                                 $(F77) -c $(FL2) -I$(INC_SRC) $< -o $@
$(OBJDIR)/otd_tools.o           :otd_tools.f;                           $(F77) -c $(FL2) -I$(INC_SRC) $< -o $@
$(OBJDIR)/linalg.o              :linalg.f;                              $(F77) -c $(FL2) -I$(INC_SRC) $< -o $@
