/*
 * This class is built for analyzing the parameters in the command line.
 * Similar to the Google FLAGS, please refer to https://code.google.com/p/gflags/
 *
 * It's easy to use. you only need to set the label, variable and default value
 * for a specific parameter. For eaxmple, if you command is as follow,
 * ./a -k 10 -evalue 0.01 -mode fast
 * then use the function of GetOption to get all the parameters. It supports int,
 * double, unsigned int, bool and string. For the above example, the functions are
 * as follow,
 *
 * Option::GetOption("-k", k, 8)
 * Option::GetOption("-evalue", evalue, 0.001)
 * Option::GetOption("mode", mode, fast)
 *
 * The first argument is the label for the particular variable, and the second argument
 * is the variable, and the third is the default value for the variable.
 *
 * In order to use class Option, you use include option.h and in the first line of your
 * main function call InitProgram(argc, argv).
 *
 * More information or suggestions please send email to haifengc at usc dot edu
 *
 */

#ifndef OPTION_H_
#define OPTION_H_

#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <string>
#include <stdio.h>

using namespace std;

class Option {
 public:
  Option();
  ~Option();

  /* argc_ and argv_ is copied from the main function. They are static variable
   *  in the class of Option so it can be used by all the funcitons in different
   *  files */
  static int argc_;
  static char** argv_;

  /* time_pre is a time label for all the output file*/
  static char* time_pre;

  /* set time_pre */
  static void GetTimePre();

  /* check a centain lable exsit or not */
  static void ChkStrExist(const char* label, bool& bVal);

  /* get the file name which is before the suffix of the file. For exmaple,
   * if the file is SRR390728.sra, then the file name is SRR390728 */
  static string GetFileName(const string& file_name);

  /* get the file extension which is after the file name. For example, if
   * the file is SRR390728.sra, then the file extension is sra */
  static string GetFileExtension(const string& file_name);

  /* get parameter for integer */
  static void GetOption(const char* label, int& intVal, const int& def);

  /* get parameter for unsigned integer */
  static void GetOption(const char* label, uint32_t & usVal,
                        const uint32_t& def);

  /* get parameter for double */
  static void GetOption(const char* label, double& doubleVal,
                        const double& def);

  /* get parameter for string, no default value, so if the command doesn't
   * indicate this parameter, the program will report the error.  */
  static void GetOption(const char* label, string& strVal);

  /* get parameter for string */
  static void GetOption(const char* label, string& strVal, const string& def);

  /* get the command which is the first argument in the command line*/
  static string GetCommand();

  /* add the time_pre to file_name */
  static string AddPreToFileName(const string& file_name);
};

/* In order to use class Option, the first line of main function should
 * be InitProgram which assigns value to argc_ and argv_ in class Option */
void InitProgram(int argc, const char* argv[]);

#endif /* OPTION_H_ */
