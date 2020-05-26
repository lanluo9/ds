disp(' ')
disp(' ')

currentEditor = matlab.desktop.editor.getActive;
originalSelection = currentEditor.Selection;

prefix = 'help NIM.';
function_name = currentEditor.SelectedText;
disp([prefix function_name])

call_help = [prefix function_name];
eval(call_help);