function modifydicomheader(input_dicom_file, output_dicom_file, f, v)

% Modify the fields of a dicom header
% 
% 2019-12-29: Created, Sam NH

% input_dicom_file = 'IM-0001-0001.dcm';
% output_dicom_file = 'IM-0001-0001-newheader.dcm';

% read data
dicom_data = dicomread(input_dicom_file);

% read header
dicom_header = dicominfo(input_dicom_file);

% write fields
for i = 1:length(f)
    dicom_header.(f{i}) = v{i};
end

% write data
dicomwrite(dicom_data, output_dicom_file, dicom_header, 'CreateMode', 'copy');

% %%
% 
% dicom_header_output = dicominfo(output_dicom_file);
% 
% dicom_header_output.SpacingBetweenSlices

