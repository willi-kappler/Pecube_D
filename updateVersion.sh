#!/bin/bash


BUILD_DATE=$(date)
BUILD_REVISION=$(git rev-parse HEAD)
OUTPUT_FILE=version.f90

BUILD_HOST=$(uname -a)
BUILD_USER=$(whoami)

echo "module m_version" > $OUTPUT_FILE
echo "  use m_logger" >> $OUTPUT_FILE
echo "  implicit none" >> $OUTPUT_FILE
echo "  contains" >> $OUTPUT_FILE
echo "    subroutine print_program_version()" >> $OUTPUT_FILE
echo "      implicit none" >> $OUTPUT_FILE
echo "" >> $OUTPUT_FILE
echo "      call log_message(\"Begin build information\")" >> $OUTPUT_FILE
echo "      call log_message(\"+- Revision: $BUILD_REVISION, User: $BUILD_USER, Date: $BUILD_DATE\")" >> $OUTPUT_FILE
echo "      call log_message(\"+- Host: $BUILD_HOST\")" >> $OUTPUT_FILE
echo "      call log_message(\"End build information\")" >> $OUTPUT_FILE
echo "" >> $OUTPUT_FILE
echo "    end subroutine print_program_version" >> $OUTPUT_FILE
echo "end module m_version" >> $OUTPUT_FILE

