import wx

version = "1"

redirect = False
redirect_file = ".libRadtran_GUI.log"

# this is not generated automatically from src_py, so I suggest not to use it.
# better refer to user guide
documentation_window = False

tooltips = True

# CE: plotting currently does not work
enable_plotting = False

# Main window
#main_window_size = 1200, 750
main_window_size = 1350, 800
#main_window_size_hints = 800, 450
main_window_size_hints = 800, 450
main_window_title = "GUI for libRadtran"

max_column_items = 20 # Maximum number of options in a single column

# Start panel
start_welcome_text = """

Welcome to the libRadtran Graphical User Interface.
        
To start a new input file, push the New button. To open an old input file, push the Open button.

Mandatory options are in red, the rest are optional.        
        
Find more information about uvspec and libRadtran at www.libradtran.org."""

# Widgets
option_style =  wx.DEFAULT_FRAME_STYLE #wx.BORDER_SIMPLE# is great for debugging!
#option_style=wx.BORDER_SIMPLE
option_background_colour = 255, -1, -1

documentation_window_position = 100, 100

plot_dpi = 96

# ColumnChoices
column_choices_y_size = 80

# ExampleList
example_thumb_size = 60, 60

maxint = 2**31 - 1 # sys.maxint, but not if python is 64 bits and wx 32 bits.

plot_button_size = 20

# Mysticthreed
threedmystic_in_GUI = True # CE this switch does not work anymore
