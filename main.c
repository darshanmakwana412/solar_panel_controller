#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sampa.h"

void calculate_irradiance_and_plot(double beta, double gamma, double albedo) {

    FILE *data_file = fopen("irradiance_data.txt", "w");
    if (data_file == NULL) {
        printf("Error opening data file for writing.\n");
        return;
    }

    sampa_data sampa;
    sampa.spa.year          = 2024;
    sampa.spa.month         = 9;
    sampa.spa.day           = 14;
    sampa.spa.timezone      = 5.5;
    sampa.spa.delta_ut1     = 0;
    sampa.spa.delta_t       = 66.4;
    sampa.spa.longitude     = 72.9151;
    sampa.spa.latitude      = 19.1331;
    sampa.spa.elevation     = 1;
    sampa.spa.pressure      = 1000;
    sampa.spa.temperature   = 30;
    sampa.spa.atmos_refract = 0.5667;
    sampa.function          = SAMPA_ALL;

    sampa.bird_sol_con = 1361.1;
    sampa.bird_ozone   = 0.3;
    sampa.bird_pwv     = 1.5;
    sampa.bird_aod     = 0.07637;
    sampa.bird_ba      = 0.85;
    sampa.bird_albedo  = albedo;

    for (int hour = 6; hour < 19; ++hour) {
        for (int minute = 0; minute < 60; minute += 1) {
            sampa.spa.hour = hour;
            sampa.spa.minute = minute;
            sampa.spa.second = 0;

            int result = sampa_calculate(&sampa);

            if (result == 0) {
                // Get the solar zenith and azimuth angles
                double zenith = sampa.spa.zenith;   // degrees
                double gamma_s = sampa.spa.azimuth; // degrees
                double DNI = sampa.dni_sul;
                double DHI = sampa.dhi_sul;

                // Convert degrees to radians
                double deg_to_rad = M_PI / 180.0;
                double alpha = 90.0 - zenith; // solar elevation angle in degrees
                double zenith_rad = zenith * deg_to_rad;
                double alpha_rad = alpha * deg_to_rad;
                double beta_rad = beta * deg_to_rad;
                double gamma_rad = gamma * deg_to_rad;
                double gamma_s_rad = gamma_s * deg_to_rad;

                // Compute cos_theta_i
                double cos_theta_i = sin(alpha_rad) * cos(beta_rad) + cos(alpha_rad) * sin(beta_rad) * cos(gamma_s_rad + gamma_rad);

                // Compute S_BEAM
                double S_BEAM = DNI * cos_theta_i;

                // Compute S_DIF
                double S_DIF = DHI * (1 + cos(beta_rad)) / 2.0;

                // Compute S_REF
                double cos_theta_z = cos(zenith_rad);
                double S_REF = albedo * (S_DIF + S_BEAM) * (1 - cos(beta_rad)) / 2.0;

                // Compute Total Irradiance
                double Total_Irradiance = fmax(0.0, S_BEAM + S_DIF + S_REF);

                // Calculate the time in hours (e.g., 13.5 for 1:30 PM)
                double time_in_hours = hour + minute / 60.0;

                // Write the time and total irradiance to the data file
                fprintf(data_file, "%.2f\t%.6f\n", time_in_hours, Total_Irradiance);
            } else {
                printf("SAMPA Error Code at %02d:%02d: %d\n", hour, minute, result);
            }
        }
    }

    fclose(data_file);
    printf("Irradiance data has been written to 'irradiance_data.txt'.\n");

    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "set terminal pngcairo size 800,600 enhanced font 'Verdana,10'\n");
        fprintf(gnuplotPipe, "set output 'irradiance_plot.png'\n");
        fprintf(gnuplotPipe, "set title 'Irradiance vs. Time'\n");
        fprintf(gnuplotPipe, "set xlabel 'Time (hours)'\n");
        fprintf(gnuplotPipe, "set ylabel 'Total Irradiance (W/m^2)'\n");
        fprintf(gnuplotPipe, "set grid\n");
        fprintf(gnuplotPipe, "plot 'irradiance_data.txt' using 1:2 with lines title 'Total Irradiance'\n");
        fflush(gnuplotPipe); // Ensure commands are sent to Gnuplot
        pclose(gnuplotPipe);
        printf("Plot has been generated and saved as 'irradiance_plot.png'.\n");
    } else {
        printf("Error: Could not open pipe to Gnuplot. Ensure Gnuplot is installed.\n");
    }
}

double compute_irradiance(double beta, double gamma, int hour, int minute) {
	double albedo = 1.0;
    sampa_data sampa;
    sampa.spa.year          = 2024;
    sampa.spa.month         = 9;
    sampa.spa.day           = 14;
    sampa.spa.timezone      = 5.5;
    sampa.spa.delta_ut1     = 0;
    sampa.spa.delta_t       = 66.4;
    sampa.spa.longitude     = 72.9151;
    sampa.spa.latitude      = 19.1331;
    sampa.spa.elevation     = 1;
    sampa.spa.pressure      = 1000;
    sampa.spa.temperature   = 30;
    sampa.spa.atmos_refract = 0.5667;
    sampa.function          = SAMPA_ALL;

    sampa.bird_sol_con = 1361.1;
    sampa.bird_ozone   = 0.3;
    sampa.bird_pwv     = 1.5;
    sampa.bird_aod     = 0.07637;
    sampa.bird_ba      = 0.85;
    sampa.bird_albedo  = 1.0;

    sampa.spa.hour = hour;
    sampa.spa.minute = minute;
    sampa.spa.second = 0;

    int result = sampa_calculate(&sampa);

    if (result == 0) {

        double zenith = sampa.spa.zenith;   // degrees
        double gamma_s = sampa.spa.azimuth; // degrees
        double DNI = sampa.dni_sul;
        double DHI = sampa.dhi_sul;


        double deg_to_rad = M_PI / 180.0;
        double alpha = 90.0 - zenith;
        double zenith_rad = zenith * deg_to_rad;
        double alpha_rad = alpha * deg_to_rad;
        double beta_rad = beta * deg_to_rad;
        double gamma_rad = gamma * deg_to_rad;
        double gamma_s_rad = gamma_s * deg_to_rad;
        double cos_theta_i = sin(alpha_rad) * cos(beta_rad) + cos(alpha_rad) * sin(beta_rad) * cos(gamma_s_rad + gamma_rad);
        double S_BEAM = DNI * cos_theta_i;
        double S_DIF = DHI * (1 + cos(beta_rad)) / 2.0;

        double cos_theta_z = cos(zenith_rad);
        double S_REF = albedo * (S_DIF + S_BEAM) * (1 - cos(beta_rad)) / 2.0;

        double Total_Irradiance = fmax(0.0, S_BEAM + S_DIF + S_REF);
        return Total_Irradiance;

    } else {
        printf("SAMPA Error Code at %02d:%02d: %d\n", hour, minute, result);
    }

}

#define DEG_TO_RAD (M_PI / 180.0)

// Define the time period and time step
#define T_START_HOUR 6
#define T_END_HOUR 18
#define TIME_STEP_MIN 5  // Time step in minutes (adjusted for computational feasibility)

// Define the beta step (Δβ)
#define DELTA_BETA 10.0  // Tilt angle step in degrees

// Define the maximum number of movements (N)
#define MAX_MOVEMENTS 5  // Adjust based on requirements

int main()
{
    double gamma = 90.0; 
    int total_minutes = (T_END_HOUR - T_START_HOUR) * 60;
    int time_steps = total_minutes / TIME_STEP_MIN + 1;

    int N_p = (int)(180.0 / DELTA_BETA) + 1;

    double *beta_positions = (double *)malloc(N_p * sizeof(double));
    for (int i = 0; i < N_p; i++)
    {
        beta_positions[i] = -90.0 + i * DELTA_BETA;
    }

    double **irradiance = (double **)malloc(N_p * sizeof(double *));
    for (int i = 0; i < N_p; i++)
    {
        irradiance[i] = (double *)malloc(time_steps * sizeof(double));
    }

    int t = 0;
    for (int hour = T_START_HOUR; hour <= T_END_HOUR; hour++)
    {
        for (int minute = 0; minute < 60; minute += TIME_STEP_MIN)
        {
            if (hour == T_END_HOUR && minute > 0)
                break;

            for (int i = 0; i < N_p; i++)
            {
                irradiance[i][t] = compute_irradiance(beta_positions[i], gamma, hour, minute);
            }
            t++;
        }
    }

    // Dynamic programming tables
    double ***dp = (double ***)malloc(time_steps * sizeof(double **));
    int ***prev_pos = (int ***)malloc(time_steps * sizeof(int **));
    for (t = 0; t < time_steps; t++)
    {
        dp[t] = (double **)malloc((MAX_MOVEMENTS + 1) * sizeof(double *));
        prev_pos[t] = (int **)malloc((MAX_MOVEMENTS + 1) * sizeof(int *));
        for (int m = 0; m <= MAX_MOVEMENTS; m++)
        {
            dp[t][m] = (double *)malloc(N_p * sizeof(double));
            prev_pos[t][m] = (int *)malloc(N_p * sizeof(int));
        }
    }

    // Initialize dp table at time 0
    for (int m = 0; m <= MAX_MOVEMENTS; m++)
    {
        for (int i = 0; i < N_p; i++)
        {
            dp[0][m][i] = irradiance[i][0];
            prev_pos[0][m][i] = -1;  // No previous position
        }
    }

    // Fill dp table
    for (t = 1; t < time_steps; t++)
    {
        for (int m = 0; m <= MAX_MOVEMENTS; m++)
        {
            for (int i = 0; i < N_p; i++)
            {
                double max_energy = -1.0;
                int best_j = -1;
                // Option 1: Stay at the same position
                dp[t][m][i] = dp[t - 1][m][i] + irradiance[i][t];
                prev_pos[t][m][i] = i;

                // Option 2: Move from another position
                if (m > 0)
                {
                    for (int j = 0; j < N_p; j++)
                    {
                        if (j != i)
                        {
                            double energy = dp[t - 1][m - 1][j] + irradiance[i][t];
                            if (energy > dp[t][m][i])
                            {
                                dp[t][m][i] = energy;
                                prev_pos[t][m][i] = j;
                            }
                        }
                    }
                }
            }
        }
    }

    double max_total_energy = -1.0;
    int best_m = -1, best_i = -1;
    for (int m = 0; m <= MAX_MOVEMENTS; m++)
    {
        for (int i = 0; i < N_p; i++)
        {
            if (dp[time_steps - 1][m][i] > max_total_energy)
            {
                max_total_energy = dp[time_steps - 1][m][i];
                best_m = m;
                best_i = i;
            }
        }
    }

    int *optimal_positions = (int *)malloc(time_steps * sizeof(int));
    int t_curr = time_steps - 1;
    int m_curr = best_m;
    int i_curr = best_i;
    while (t_curr >= 0)
    {
        optimal_positions[t_curr] = i_curr;
        int i_prev = prev_pos[t_curr][m_curr][i_curr];
        if (i_prev != i_curr)
        {
            m_curr--;
        }
        i_curr = i_prev;
        t_curr--;
    }


    printf("Time\tBeta\n");
    t = 0;
    int prev_position = -1;
    for (int hour = T_START_HOUR; hour <= T_END_HOUR; hour++)
    {
        for (int minute = 0; minute < 60; minute += TIME_STEP_MIN)
        {
            if (hour == T_END_HOUR && minute > 0)
                break;

            int position = optimal_positions[t];
            if (position != prev_position)
            {
                printf("%02d:%02d\t%.1f\n", hour, minute, beta_positions[position]);
                prev_position = position;
            }
            t++;
        }
    }

    printf("Total Energy Received: %.2f Wh/m^2\n", max_total_energy * TIME_STEP_MIN / 60.0);
    for (t = 0; t < time_steps; t++)
    {
        for (int m = 0; m <= MAX_MOVEMENTS; m++)
        {
            free(dp[t][m]);
            free(prev_pos[t][m]);
        }
        free(dp[t]);
        free(prev_pos[t]);
    }
    free(dp);
    free(prev_pos);

    for (int i = 0; i < N_p; i++)
    {
        free(irradiance[i]);
    }
    free(irradiance);
    free(beta_positions);
    free(optimal_positions);

    return 0;
}
