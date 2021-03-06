#
# This file is part of GRAMPC - (https://sourceforge.net/projects/grampc/)
#
# GRAMPC -- A software framework for embedded nonlinear model predictive
# control using a gradient-based augmented Lagrangian approach
#
# Copyright (C) 2014-2018 by Andreas Voelz, Bartosz Kaepernick,
# Felix Mesmer, Knut Graichen, Soenke Rhein, Tobias Englert. Developed at the 
# Institute of Measurement, Control, and Microtechnology, Ulm University.
# All rights reserved.
#
# GRAMPC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as 
# published by the Free Software Foundation, either version 3 of 
# the License, or (at your option) any later version.
#
# GRAMPC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public 
# License along with GRAMPC. If not, see <http://www.gnu.org/licenses/>.
#


#
# compiler
#
COMPILER = gcc
AR       = ar r
RM       = rm -f
MKDIR    = mkdir -p
OBJEXT   = o
LIBEXT   = a


#
# path of header and source files
#
HEADER_PATH = ./include
SOURCE_PATH = ./src
LIBS_PATH   = ./libs


#
# flags and settings
#
HEADER = -I${HEADER_PATH}
CFLAGS = -c -O3 -Wall
LIBS   = 


#
# object files and lib
#		
GRAMPC_OBJS = \
		$(SOURCE_PATH)/euler1.$(OBJEXT) \
		$(SOURCE_PATH)/eulermod2.$(OBJEXT) \
		$(SOURCE_PATH)/heun2.$(OBJEXT) \
		$(SOURCE_PATH)/ruku45.$(OBJEXT) \
		$(SOURCE_PATH)/rodas.$(OBJEXT) \
        $(SOURCE_PATH)/trapezodial.$(OBJEXT) \
        $(SOURCE_PATH)/simpson.$(OBJEXT) \
		$(SOURCE_PATH)/grampc_init.$(OBJEXT) \
		$(SOURCE_PATH)/grampc_mess.$(OBJEXT) \
		$(SOURCE_PATH)/grampc_run.$(OBJEXT) \
		$(SOURCE_PATH)/grampc_setopt.$(OBJEXT) \
		$(SOURCE_PATH)/grampc_setparam.$(OBJEXT) \
		$(SOURCE_PATH)/grampc_util.$(OBJEXT) 
GRAMPC_LIB = \
		$(LIBS_PATH)/libgrampc.$(LIBEXT)


#
# targets
#
all: $(GRAMPC_LIB)

$(GRAMPC_LIB): $(GRAMPC_OBJS)
			@echo 'Building library: $@'
			$(MKDIR) $(LIBS_PATH)
			$(AR) -o "$@" $(GRAMPC_OBJS)
			@echo 'Finished building: $@'
			@echo ''

$(GRAMPC_OBJS): %.$(OBJEXT): %.c
			@echo 'Building file: $@'
			$(COMPILER) -o "$@" "$<" $(CFLAGS) $(HEADER) $(LIBS) 
			@echo 'Finished building: $@'
			@echo ''

clean: 
			@echo 'Cleaning up built files:'
			$(RM) $(GRAMPC_LIB) $(GRAMPC_OBJS)
			@echo ''

#
# end of make
#		
