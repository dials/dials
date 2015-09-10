import subprocess
import wx


def create_gui_dat():
  class element(object):
    def __init__(self):
      self.lbl = "\n Label here \n"
      self.cmd = "echo \"a default command line goes here\""

  local_lst = []

  import_btn = element()
  import_btn.lbl = "\n import \n"
  import_btn.cmd = "dials.import"
  local_lst.append(import_btn)

  find_btn = element()
  find_btn.lbl = "\n find_spots \n"
  find_btn.cmd = "dials.find_spots"
  local_lst.append(find_btn)

  index_btn = element()
  index_btn.lbl = "\n index \n"
  index_btn.cmd = "dials.index"
  local_lst.append(index_btn)

  dummy_tst = element()
  local_lst.append(dummy_tst)

  return local_lst


class MyWidg(wx.Frame):
  def __init__(self, parent, id, title, data_in):
    super(MyWidg, self).__init__(parent, id, title)

    self.step_btn_lst = data_in

    for step_btn in self.step_btn_lst:
      exmpl_btt = wx.Button(self, -1, step_btn.lbl, (0,0))
      step_btn.wx_btn = exmpl_btt
      step_btn.bt_id = exmpl_btt.GetId()

    main_panel = wx.Panel(self)

    go_btn = wx.Button(self, -1, " Go ", (0,0))
    self.cli_txt = wx.TextCtrl(self, -1, size = (800,20), style = wx.DEFAULT)


    vbox = wx.BoxSizer(wx.VERTICAL)

    for step_btn in self.step_btn_lst:
      vbox.Add(step_btn.wx_btn, 0, border = 3)


    hbox = wx.BoxSizer(wx.HORIZONTAL)
    hbox.Add(vbox, 0, border = 3)
    hbox.Add(main_panel, 0, border = 3)

    bt_box = wx.BoxSizer(wx.HORIZONTAL)
    bt_box.Add(self.cli_txt,wx.EXPAND | wx.ALL, border = 3)
    bt_box.Add(go_btn, 0, border = 3)

    bg_box = wx.BoxSizer(wx.VERTICAL)
    bg_box.Add(hbox, 0, border = 3)
    bg_box.Add(bt_box, 0, border = 3)

    self.SetSizer(bg_box)
    bg_box.Fit(self)

    for step_btn in self.step_btn_lst:
      step_btn.wx_btn.Bind(wx.EVT_BUTTON, self.on_btn)

    go_btn.Bind(wx.EVT_BUTTON, self.on_go_button)


  def on_btn(self, event):

    for step_btn in self.step_btn_lst:
      if(step_btn.bt_id == event.GetId()):
        self.cli_txt.SetValue(step_btn.cmd)


  def on_go_button(self, event):
    txt_to_run = self.cli_txt.GetValue()
    subprocess.call(txt_to_run, shell=True)


class MyApp(wx.App):
  def OnInit(self):

    gui_dat = create_gui_dat()

    frame = MyWidg(parent = None, id = -1, title = 'Main Frame', data_in = gui_dat)
    frame.Show(True)
    return True

app = MyApp(0)
app.MainLoop()
