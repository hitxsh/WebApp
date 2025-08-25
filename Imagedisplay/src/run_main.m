function run_main(inputPath, outputPath)
    % Copy input to where your script expects it
    copyfile(inputPath, '/Users/aayush_kumar/Documents/PhD/Experiments/May 2022 onwards/Shape Characterization Codes/temp_input.tif');

    % Run your original script
    Codemsir;

    % Assume output is saved to a known name, e.g., temp_output.tif
    % (Adjust this part if your output filename is different)
    copyfile('/Users/aayush_kumar/Documents/PhD/Experiments/May 2022 onwards/Shape Characterization Codes/output.tif', outputPath);
end
