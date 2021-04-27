# EEGBrowser
The EEGBrowser is intended to be a modern drop-in replacement for EEGLAB's `eegplot` function based on [MoBILAB's](https://sccn.ucsd.edu/wiki/MoBILAB) visualization functions.

![EEGBrowser](https://github.com/aojeda/EEGBrowser/blob/master/resources/snapshot.png)

## Install
* [Download](https://github.com/aojeda/EEGBrowser/archive/master.zip)
* Uncompress
* Rename the folder as "EEGBrowser", and place it in your `eeglab/plugins` folder.
* Run `eeglab` or `eeglab redraw` (if you don't want to clean your current EEG structure)

![Plugin](https://github.com/aojeda/EEGBrowser/blob/master/resources/snapshot2.png)

**Note:** Although EEG Browser is based on [MoBILAB](https://sccn.ucsd.edu/wiki/MoBILAB), the latter is not needed for it to work.

## Example

Plot channel data (default)
```matlab
pop_eegbrowser(EEG);
```
or
```matlab
pop_eegbrowser(EEG,1);
```
Plot IC activations:
```matlab
pop_eegbrowser(EEG,0);
```

To get the handle to the widget (useful for adding your own customizations):
```matlab
hBrowser = pop_eegbrowser(EEG,0);
```

Other features:
* Use the `<<` and `>>` keys to move to the next page centered around the selected event marker.
* Use `-` or `+` buttons to reduce or increase the scale respectively.
* Click on the figure and use `-` or `+` keys in your keypad to reduce or increase the scale respectively.
* Use the mouse to select bad trials. See example [here](https://github.com/aojeda/EEGBrowser/wiki/Example:-manual-trial-rejection).

## Developers

If you want to extend the capabilities of the EEG Browser follow these steps:
* Create a class inheriting from `EEGBrowser`, let's call it `MyEEGBrowser`. See example below.
* Then create your own `pop_myeegbrowser` or modify `pop_eegbrowser` [here](https://github.com/aojeda/EEGBrowser/blob/master/pop_eegbrowser.m#L30) to call your class passing in the EEG structure.

### Example of creating a new browser class ###
``` matlab
classdef MyEEGBrowser < EEGBrowser
    properties
      % Declare the properties that you need here for plotting.
      % Properties store data that is used by the class for plotting.
    end

    % Constructor
    function obj = MyEEGBrowser(EEG)

      % Call EEGBrowser's constructor
      obj@EEGBrowser(EEG);

      %--
      % Do your thing here. You can for instance add buttons, change default
      % settings, and so on.
      hButton = uicontrol('Parent', obj.figureHandle,...
                          'Style', 'pushbutton',...
                          'TooltipString','This does something',...
                          'Units','Normalized',...
                          'Callback',@onPush);
      %--
    end
end

% Define your callbacks
function onPush(src,evnt)
  obj = src.Parent.UserData;  % Retrieve your MyEEGBrowser object
                              % so that you can modify its properties
  %--
  % Implement the callback here.
  %--
end
```

Save the code above in *MyEEGBrowser.m* and that's it. See [this](https://www.mathworks.com/help/matlab/matlab_oop/create-a-simple-class.html) example for learning more about MATLAB classes.

## Version history
* v1.1 - fix issue associated with modified EEGLAB menu in development branch
* v1.0 - initial release

