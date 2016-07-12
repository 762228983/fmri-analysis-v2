function unix_freesurfer_version(freesurfer_version, command)

global root_directory;

freesurfer_directory = [root_directory '/freesurfer'];
setupstr = [...
    'export FREESURFER_HOME=/software/Freesurfer/' freesurfer_version ';'  ...
    ' source ${FREESURFER_HOME}/SetUpFreeSurfer.sh;'...
    ' export SUBJECTS_DIR=' freesurfer_directory ';'];

myfsaverage_directory = [freesurfer_directory '/myfsaverage'];
if ~exist(myfsaverage_directory, 'dir')
    mkdir(myfsaverage_directory);
    unix(['cp -r /mindhive/nklab/u/svnh/freesurfer/myfsaverage/surf ' myfsaverage_directory '/surf']);
    unix(['cp -r /mindhive/nklab/u/svnh/freesurfer/myfsaverage/label ' myfsaverage_directory '/label']);
    unix(['cp -r /mindhive/nklab/u/svnh/freesurfer/myfsaverage/mri ' myfsaverage_directory '/mri']);
    unix(['cp -r /mindhive/nklab/u/svnh/freesurfer/myfsaverage/mri ' myfsaverage_directory '/mri.2mm']);
    unix(['cp -r /mindhive/nklab/u/svnh/freesurfer/myfsaverage/scripts ' myfsaverage_directory '/scripts']);
end

unix([setupstr command]);
