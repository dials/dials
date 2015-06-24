import subprocess
import wx

class MyWidg(wx.Frame):
    def __init__(self, parent, id, title):
        super(MyWidg, self).__init__(parent, id, title)
        statusbar = self.CreateStatusBar()
        statusbar.SetStatusText('Running')

        frs_btt = wx.Button(self, -1, "\n import \n", (0,0))
        seq_btt = wx.Button(self, -1, "\n find_spots \n", (0,0))
        trd_btt = wx.Button(self, -1, "\n index \n", (0,0))
        main_panel = wx.Panel(self)
        self.cli_txt = wx.TextCtrl(self, -1, size=(1500,20), style=wx.DEFAULT)

        go_btt = wx.Button(self, -1, " Go ", (0,0))

        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(frs_btt, 0, border = 3)
        vbox.Add(seq_btt, 0, border = 3)
        vbox.Add(trd_btt, 0, border = 3)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(vbox, 0, border = 3)
        hbox.Add(main_panel, 0, border = 3)

        bt_box = wx.BoxSizer(wx.HORIZONTAL)
        bt_box.Add(self.cli_txt,wx.EXPAND | wx.ALL, border = 3)
        bt_box.Add(go_btt, 0, border = 3)

        bg_box = wx.BoxSizer(wx.VERTICAL)
        bg_box.Add(hbox, 0, border = 3)
        bg_box.Add(bt_box, 0, border = 3)

        self.SetSizer(bg_box)
        bg_box.Fit(self)

        frs_btt.Bind(wx.EVT_BUTTON, self.on_imp_button)
        seq_btt.Bind(wx.EVT_BUTTON, self.on_fnd_button)
        trd_btt.Bind(wx.EVT_BUTTON, self.on_ind_button)
        go_btt.Bind(wx.EVT_BUTTON, self.on_go_button)

    def on_imp_button(self, event):
        self.cli_txt.SetValue("dials.import")

    def on_fnd_button(self, event):
        self.cli_txt.SetValue("dials.find_spots")

    def on_ind_button(self, event):
        self.cli_txt.SetValue("dials.index")

    def on_go_button(self, event):
        txt_to_run = self.cli_txt.GetValue()
        #'print "string to run =", txt_to_run
        subprocess.call(txt_to_run, shell=True)


class MyApp(wx.App):
    def OnInit(self):
        frame = MyWidg(None, -1, 'menu1.py')
        frame.Show(True)
        return True

app = MyApp(0)
app.MainLoop()
