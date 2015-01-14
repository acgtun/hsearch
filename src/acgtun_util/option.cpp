/* The detail description of each function please refer to head file */

#include "option.hpp"
#include <time.h>

/* argc_ and argv_ are static class members, they should be initiated
 * when the program starts in order to assign memory for these two variable  */
int Option::argc_ = 0;
char** Option::argv_ = NULL;
char* Option::time_pre = NULL;

Option::Option() {
}
Option::~Option() {
}

void Option::ChkStrExist(const char* label, bool& bVal) {
  for (int i = 1; i < argc_; i++) {
    if (strcmp(argv_[i], label) == 0) {
      bVal = true;
      return ;
    }
  }
  bVal = false;
}

string Option::GetCommand() {
  string strVal = argv_[0];
  int pos = strVal.find_last_of("/");
  string command = strVal.substr(pos + 1);
  return command;
}

string Option::GetFileName(const string& file_name) {
  uint32_t p = file_name.find_last_of('.');

  return file_name.substr(0, p);
}

string Option::GetFileExtension(const string& file_name) {
  printf("%s\n", file_name.c_str());
  uint32_t p = file_name.find_last_of('.');

  return file_name.substr(p);
}

void Option::GetOption(const char* label, int& intVal, const int& def) {
  for (int i = 1; i < argc_; i++) {
    if (strcmp(argv_[i], label) == 0 && i + 1 <= argc_ - 1) {
      if (argv_[i + 1][0] != '-') {
        intVal = atoi(argv_[i + 1]);
        return;
      }
    }
  }
  intVal = def;
}

void Option::GetOption(const char* label, uint32_t & usVal,
                       const uint32_t& def) {
  for (int i = 1; i < argc_; i++) {
    if (strcmp(argv_[i], label) == 0 && i + 1 <= argc_ - 1) {
      if (argv_[i + 1][0] != '-') {
        usVal = atoi(argv_[i + 1]);
        return;
      }
    }
  }
  usVal = def;
}

void Option::GetOption(const char* label, double& doubleVal,
                       const double& def) {
  for (int i = 1; i < argc_; i++) {
    if (strcmp(argv_[i], label) == 0 && i + 1 <= argc_ - 1) {
      if (argv_[i + 1][0] != '-') {
        doubleVal = strtod(argv_[i + 1], NULL);
        return;
      }
    }
  }
  doubleVal = def;
}

void Option::GetOption(const char* label, string& strVal) {
  for (int i = 1; i < argc_; i++) {
    if (strcmp(argv_[i], label) == 0 && i + 1 <= argc_ - 1) {
      if (argv_[i + 1][0] != '-') {
        strVal = argv_[i + 1];
        return;
      }
    }
  }
  printf("please set parameter %s", label);
  exit (EXIT_FAILURE);
}

void Option::GetOption(const char* label, string& strVal, const string& def) {
  for (int i = 1; i < argc_; i++) {
    if (strcmp(argv_[i], label) == 0 && i + 1 <= argc_ - 1) {
      if (argv_[i + 1][0] != '-') {
        strVal = argv_[i + 1];
        return;
      }
    }
  }
  strVal = def;
}

string Option::AddPreToFileName(const string& file_name) {
  string ret;
  ret = Option::time_pre;
  ret += file_name;

  return ret;
}

void Option::GetTimePre() {
  time_t rawtime;
  struct tm * timeinfo;

  time(&rawtime);
  timeinfo = localtime(&rawtime);

  char pre[100] = { '\0' };
  sprintf(pre, "%04d-%02d-%02d-%02d-%02d-%02d_", 1900 + timeinfo->tm_year,
          1 + timeinfo->tm_mon, timeinfo->tm_mday, timeinfo->tm_hour,
          timeinfo->tm_min, timeinfo->tm_sec);
  Option::time_pre = (char *) malloc(sizeof(char) * (strlen(pre) + 1));
  strcpy(time_pre, pre);
}

void InitProgram(int argc, const char* argv[]) {
  Option::argc_ = argc;
  Option::argv_ = (char**) malloc(sizeof(char *) * argc);
  for (int i = 0; i < argc; i++) {
    Option::argv_[i] = (char *) malloc(sizeof(char) * (strlen(argv[i]) + 1));
    strcpy(Option::argv_[i], argv[i]);
  }
  Option::GetTimePre();
}
