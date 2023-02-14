/* 4D检验法去除异常数据 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
// #define MIN_TIME 0.000001
// #define DEBUG 1

int main(int argc, char *argv[]) {
  int flag = 0, count = 0, row = 0, base = 0, i = 0, j = 0;
  FILE *data1;
  long double time_base[100];
  long double time_base_sum = 0;
  int start, end;
  long double time_sum = 0;
  long double time_output_sum = 0;
  int time_output_count = 0;

#ifdef DEBUG
  printf("4D Test will opening file(%s).\r\n", argv[1]);
#endif // DEBUG
  data1 = fopen(argv[1], "r");
  if (!data1) {
    printf("Error in opening file(%s).\r\n", argv[1]);
    exit(1);
  }

  row = 0;
  while (fscanf(data1, "%Lf\n", &time_base[row]) != EOF)
    row++;
#ifdef DEBUG
  printf("the file(%s) total has %d lines data.\r\n", argv[1], row);
#endif // DEBUG

  start = 0;
  end = row - start;
  for (i = start; i < end; i++) {
    time_sum += time_base[i];
  }
#ifdef DEBUG
  printf("the sum of the %d lines data is %.10LF.\r\n", row, time_sum);
#endif // DEBUG

  for (i = start; i < end; i++) {
#ifdef DEBUG
    printf("the %dth data is %.10LF.\r\n", i, time_base[i]);
#endif // DEBUG

    /*
    // 设置时间最小值，如果小于该值直接舍弃.
    if (time_base[i] < MIN_TIME) {
      continue;
    }*/
    time_base_sum = time_sum - time_base[i];
#ifdef DEBUG
    printf("time_base_sum is %.10LF.\r\n", time_base_sum);
#endif // DEBUG

    long double avg1 = time_base_sum / (row - 1);
#ifdef DEBUG
    printf("avg1 is %.10LF.\r\n", avg1);
#endif // DEBUG
    
    long double avg2 = 0.0;
    for (j = start; j < end; j++) {

      if (j != i) {   
        avg2 += (fabsl(time_base[j] - avg1));
      }
    }
#ifdef DEBUG
    printf("avg2 is %.10LF.\r\n", avg2);
#endif // DEBUG

    avg2 = (avg2 / (row - 1)) * 4;
    if (0 == (time_base[i] - avg1) && 0 == avg2) {
      time_output_sum += time_base[i];
      time_output_count++;
      continue;
    }
    if (fabsl(time_base[i] - avg1) < avg2) {
      time_output_sum += time_base[i];
      time_output_count++;
    }
#ifdef DEBUG
    printf("after calculate the %dth data(%.10LF),the count is %d .\r\n", i,
           time_base[i], time_output_count);
#endif // DEBUG
  }

  // printf("%s\t", argv[1]);
  if (time_output_count != 0) {
    printf("%.10Lf\n", time_output_sum / time_output_count);
#ifdef DEBUG
    printf("All the count is %d .\r\n", time_output_count);
#endif // DEBUG
  } else {
    printf("Error:input numbers are all incorrect!\n");
  }
  fclose(data1);

  return (0);
}
