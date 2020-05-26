%% Mac version
disp('')
disp('')
currentEditor = matlab.desktop.editor.getActive; 
originalSelection = currentEditor.Selection; 
assert(originalSelection(1)==originalSelection(3)); %Check that multiple lines are not selected 
currentEditor.Selection = [originalSelection(1) 1 originalSelection(1) Inf]; %Select the whole line 
eval(currentEditor.SelectedText); %Run the whole line 
currentEditor.Selection = originalSelection; %Reset selection to original state 
clear currentEditor originalselection

%% windows version seems to be the same
currentEditor = matlab.desktop.editor.getActive;
originalSelection = currentEditor.Selection;
assert(originalSelection(1)==originalSelection(3));%Check that multiple lines are not selected
currentEditor.Selection = [originalSelection(1) 1 originalSelection(1) Inf];%Select the whole line
eval(currentEditor.SelectedText);%Run the whole line
currentEditor.Selection = originalSelection;%Reset selection to original state
clear currentEditor originalselection

%% windows version, improved to go to next line automatically 
currentEditor = matlab.desktop.editor.getActive;
originalSelection = currentEditor.Selection;
assert(originalSelection(1)==originalSelection(3));%Check that multiple lines are not selected
currentEditor.Selection = [originalSelection(1) 1 originalSelection(1) Inf];%Select the whole line
disp(currentEditor.SelectedText)
eval(currentEditor.SelectedText);%Run the whole line
currentEditor.Selection = originalSelection + 1;%Reset selection to original state + 1 line (go to next line)
clear currentEditor originalselection

%% fix shifting bug
disp('')
disp('')
currentEditor = matlab.desktop.editor.getActive;
originalSelection = currentEditor.Selection;
assert(originalSelection(1)==originalSelection(3));%Check that multiple lines are not selected
currentEditor.Selection = [originalSelection(1) 1 originalSelection(1) Inf];%Select the whole line
disp(currentEditor.SelectedText)
eval(currentEditor.SelectedText);%Run the whole line
% currentEditor.Selection = originalSelection + 1;%Reset selection to original state + 1 line (go to next line)
currentEditor.Selection = [originalSelection(1) + 1, 1, originalSelection(1) + 1, 1];%Reset selection to original state + 1 line (go to next line)
clear currentEditor originalselection