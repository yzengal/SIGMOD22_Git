#include "monitor.h"

void save_time(program_t& prog) {
	// clock_gettime(CLOCK_MONOTONIC, &prog.real);
	// clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &prog.proc);
	gettimeofday(&prog.real, NULL);
}

double calc_time(const program_t& st, const program_t& ed) {
	double ret;

	// #ifdef USING_PROC_TIME
	// double proc_st = st.proc.tv_sec + (float)st.proc.tv_nsec / 1e9;
	// double proc_ed = ed.proc.tv_sec + (float)ed.proc.tv_nsec / 1e9;
	// ret = proc_ed - proc_st;
	// #else
	// double real_st = st.real.tv_sec + (float)st.real.tv_nsec / 1e9;
	// double real_ed = ed.real.tv_sec + (float)ed.real.tv_nsec / 1e9;
	ret = (ed.real.tv_sec - st.real.tv_sec) + (ed.real.tv_usec - st.real.tv_usec) / 1000000.0;
	// #endif 

	return ret;
}

int get_proc_status(int pid, const char * mark) {
	const int BUFFER_SIZE = 5120;
	FILE * pf;
	char fn[BUFFER_SIZE], buf[BUFFER_SIZE];
	int ret = 0;

	sprintf(fn, "/proc/%d/status", pid);
	pf = fopen(fn, "re");
	int m = strlen(mark);
	while (pf && fgets(buf, BUFFER_SIZE - 1, pf)) {

		buf[strlen(buf) - 1] = 0;
		if (strncmp(buf, mark, m) == 0) {
			sscanf(buf + m + 1, "%d", &ret);
		}
	}

	if (pf)
		fclose(pf);

	return ret;
}

// usedMemory: KByte, usedTime: second
void watchSolution(pid_t pid, int& usedMemory, int& usedTime, int limitTime) {
	int tempMemory;
	int status;
	struct rusage ruse;

	usedMemory= get_proc_status(pid, "VmRSS:");

	while (1) {
		// check the usage
		wait4(pid, &status, 0, &ruse);

		tempMemory = get_proc_status(pid, "VmPeak:");
		if (tempMemory > usedMemory)
			usedMemory = tempMemory;

		if (WIFEXITED(status))
			break;

		if (WIFSIGNALED(status)) {
			printf("program is terminated unexpectedly.\n");
			break;
		}
	}

	usedTime += (ruse.ru_utime.tv_sec + ruse.ru_utime.tv_usec / 1000000);
	usedTime += (ruse.ru_stime.tv_sec + ruse.ru_stime.tv_usec / 1000000);
}

void watchSolutionOnce(pid_t pid, int& usedMemory) {
	int tempMemory;
	int status;
	struct rusage ruse;

	tempMemory = get_proc_status(pid, "VmRSS:");
	if (tempMemory > usedMemory)
		usedMemory = tempMemory;
	
	wait4(pid, &status, 0, &ruse);
	tempMemory = get_proc_status(pid, "VmPeak:");
	if (tempMemory > usedMemory)
		usedMemory = tempMemory;
}