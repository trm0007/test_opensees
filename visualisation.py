import pandas as pd
from import_ import *
from units import *

def visualization():
  # Visualization
  # Pushover Curve
  df_R_X = pd.read_table('FGU_RC3DF_files/Pushover_Horizontal_ReactionsX.out', sep = " ", header = None,
                        names=["Pseudo-Time","R1_X","R2_X","R3_X","R4_X"])
  df_R_Y = pd.read_table('FGU_RC3DF_files/Pushover_Horizontal_ReactionsY.out', sep = " ", header = None,
                        names=["Pseudo-Time","R1_Y","R2_Y","R3_Y","R4_Y"])

  df_R_X['sum_R'] = df_R_X.values[:,1:5].sum(axis =1)
  df_R_Y['sum_R'] = df_R_Y.values[:,1:5].sum(axis =1)

  df_D_X = pd.read_table('FGU_RC3DF_files/Pushover_Story_DisplacementX.out', sep = " ", header = None,
                        names=["Pseudo-Time","D1_X","D2_X","D3_X","D4_X"])
  df_D_Y = pd.read_table('FGU_RC3DF_files/Pushover_Story_DisplacementY.out', sep = " ", header = None,
                        names=["Pseudo-Time","D1_Y","D2_Y","D3_Y","D4_Y"])

  df_D_X['avg_D'] = df_D_X.values[:,1:5].mean(axis = 1)
  df_D_Y['avg_D'] = df_D_Y.values[:,1:5].mean(axis = 1)

  plt.figure(figsize=(10,5))

  plt.plot(df_D_X['avg_D'], -df_R_X['sum_R'], color = '#C0392B', linewidth=1.5)
  plt.plot(df_D_Y['avg_D'], -df_R_Y['sum_R'], color = '#27AE60', linewidth=1.5)


  plt.ylabel('Base Shear (KN)', {'fontname':'Cambria', 'fontstyle':'italic','size':14})
  plt.xlabel('Average of Roof Displacement (m)', {'fontname':'Cambria', 'fontstyle':'italic','size':14})
  plt.grid(which='both')
  plt.title('Pushover Curve',{'fontname':'Cambria', 'fontstyle':'normal','size':16})
  plt.yticks(fontname = 'Cambria', fontsize = 14)
  plt.xticks(fontname = 'Cambria', fontsize = 14);
  plt.legend(['X-Direction', 'Y-Direction'],prop={'family':'Cambria','size':14});

  # Ground Motion histroy
  G_M =np.loadtxt('FGU_RC3DF_files/acc_1.txt')
  times = np.arange(0,0.02*len(G_M),0.02)
  plt.figure(figsize=(12,4))
  plt.plot(times,G_M, color = '#6495ED', linewidth=1.2)
  plt.ylabel('Acceleration (m/s2)', {'fontname':'Cambria', 'fontstyle':'italic','size':14})
  plt.xlabel('Time (sec)', {'fontname':'Cambria', 'fontstyle':'italic','size':14})
  plt.grid(which='both')
  plt.title('Time history of Ground Motion record',{'fontname':'Cambria', 'fontstyle':'normal','size':16})
  plt.yticks(fontname = 'Cambria', fontsize = 14);

  # Time history of displacement and acceleration
  story_disp_X = np.loadtxt('FGU_RC3DF_files/TimeHistory_Story_DisplacementX1.1.out')
  story_disp_Y = np.loadtxt('FGU_RC3DF_files/TimeHistory_Story_DisplacementY1.1.out')

  plt.figure(figsize=(12,5))
  plt.plot(story_disp_X[:,0], story_disp_X[:,1], color = '#DE3163', linewidth=1.2)
  plt.plot(story_disp_Y[:,0], story_disp_Y[:,2], color = '#FFBF00', linewidth=1.2)
  plt.ylabel('Horizontal Displacement (m)', {'fontname':'Cambria', 'fontstyle':'italic','size':14})
  plt.xlabel('Time (sec)', {'fontname':'Cambria', 'fontstyle':'italic','size':14})
  plt.grid(which='both')
  plt.title('Time history of horizontal dispacement',{'fontname':'Cambria', 'fontstyle':'normal','size':16})
  plt.yticks(fontname = 'Cambria', fontsize = 14);
  plt.xticks(fontname = 'Cambria', fontsize = 14);
  plt.legend(['X-Direction', 'Y-Direction'], prop={'family':'Cambria','size':14});

  story_accel_X = np.loadtxt('FGU_RC3DF_files/TimeHistory_Story_AccelerationX1.1.out')
  story_accel_Y = np.loadtxt('FGU_RC3DF_files/TimeHistory_Story_AccelerationY1.1.out')

  plt.figure(figsize=(12,5))
  plt.plot(story_accel_X[:,0], story_accel_X[:,1], color = '#DE3163', linewidth=1.2)
  plt.plot(story_accel_Y[:,0], story_accel_Y[:,2], color = '#FFBF00', linewidth=1.2)
  plt.ylabel('Horizontal Acceleration (m/s2)', {'fontname':'Cambria', 'fontstyle':'italic','size':14})
  plt.xlabel('Time (sec)', {'fontname':'Cambria', 'fontstyle':'italic','size':14})
  plt.grid(which='both')
  plt.title('Time history of horizontal acceleration',{'fontname':'Cambria', 'fontstyle':'normal','size':16})
  plt.yticks(fontname = 'Cambria', fontsize = 14);
  plt.xticks(fontname = 'Cambria', fontsize = 14);
  plt.legend(['X-Direction', 'Y-Direction'], prop={'family':'Cambria','size':14});
  # plt.show()
