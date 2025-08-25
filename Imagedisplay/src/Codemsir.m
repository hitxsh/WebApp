function Codemsir(filename)
    save('input_temp.mat', 'filename');  % Save argument to file
    Codemsir_internal();                 % Call the real processing script
end
