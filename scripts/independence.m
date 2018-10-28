% MATLAB script to check independence of experimental effects.
% If they are not linearly independent, they cannot be used
% together in a multiple regression analysis.

% Load experimental effects.
n=input('Number of effects of interest? ');
disp('')
disp('')
disp('******************************************************')

X=[];
for i=1:n
    disp('Type in the name for the textfile')
    disp(['for experimental effect #' num2str(i) ','])
    disp('WITHOUT .txt extension.')
    filename=input('Textfile name: ','s');
    eval(['load ' filename '.txt'])
    eval([filename '=' filename '-mean(' filename ');']);
    eval(['X=[X ' filename '];'])
    disp('')
    disp('')
    disp('******************************************************')
end

% Test invertibility.
disp('Here is the determinant of X-transpose times X:')
disp('')
disp(det(X'*X))
disp('If it is zero, then your independent variables are')
disp('NOT linearly independent, and you cannot use them')
disp('together in a multiple regression analysis.')
