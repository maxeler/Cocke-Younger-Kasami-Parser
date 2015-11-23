CPU_CODE_DIR="../../CPUCode"

# delegates any target to the CPU Code Makefile
%:
	make -C $(CPU_CODE_DIR) RUNRULE="Vectis" $@

# default target
all:
	make -C $(CPU_CODE_DIR) RUNRULE="Vectis"

installer_no_dependencies:
	echo "This target cannot be handled by this version of the makefile system" > /dev/stderr
