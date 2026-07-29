import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from twixtools import read_twix, map_twix

def phyllotaxis3D(m_lNumberOfFrames, m_lProjectionsPerFrame, flagSelf):
    NProj = m_lNumberOfFrames * m_lProjectionsPerFrame
    lTotalNumberOfProjections = NProj  # For UTE, where we are going from pole to pole of the sphere, lTotalNumberOfProjections = NProj

    m_adAzimuthalAngle = np.zeros(NProj)
    m_adPolarAngle = np.zeros(NProj)
    x = np.zeros(NProj)
    y = np.zeros(NProj)
    z = np.zeros(NProj)

    if flagSelf:
        kost = np.pi / (2 * np.sqrt(lTotalNumberOfProjections - m_lNumberOfFrames))
    else:
        kost = np.pi / (2 * np.sqrt(lTotalNumberOfProjections))

    Gn = (1 + np.sqrt(5)) / 2
    Gn_ang = 2 * np.pi - (2 * np.pi / Gn)
    # Gn_ang = (2*pi / Gn)
    count = 1

    for lk in range(1, m_lProjectionsPerFrame + 1):
        for lFrame in range(1, m_lNumberOfFrames + 1):

            linter = lk + (lFrame - 1) * m_lProjectionsPerFrame

            if flagSelf and lk == 1:
                m_adPolarAngle[linter - 1] = 0
                m_adAzimuthalAngle[linter - 1] = 0
            else:
                m_adPolarAngle[linter - 1] = kost * np.sqrt(count)
                m_adAzimuthalAngle[linter - 1] = (count) * Gn_ang % (2 * np.pi)
                count = count + 1

            x[linter - 1] = np.sin(m_adPolarAngle[linter - 1]) * np.cos(m_adAzimuthalAngle[linter - 1])
            y[linter - 1] = np.sin(m_adPolarAngle[linter - 1]) * np.sin(m_adAzimuthalAngle[linter - 1])
            z[linter - 1] = np.cos(m_adPolarAngle[linter - 1])

    return m_adPolarAngle, m_adAzimuthalAngle, x, y, z

def computePhyllotaxis(N, nseg, nshot, flagSelfNav, flagPlot, flagRosettaTraj=False):
    def sph2cart(az, elev, r):
        """
        Transform spherical to Cartesian coordinates.

        Parameters:
        - az: Azimuth angle (in radians)
        - elev: Elevation angle (in radians)
        - r: Radius

        Returns:
        - x: x-coordinate in Cartesian space
        - y: y-coordinate in Cartesian space
        - z: z-coordinate in Cartesian space
        """
        z = r * np.sin(elev)
        rcoselev = r * np.cos(elev)
        x = rcoselev * np.cos(az)
        y = rcoselev * np.sin(az)
        return x, y, z
    
    polarAngle, azimuthalAngle, vx, vy, vz = phyllotaxis3D(nshot, nseg, flagSelfNav)

    r = np.linspace(-0.5, 0.5 - (1/N), N)
    azimuthal = np.tile(azimuthalAngle, (N, 1))
    polar = np.tile(np.pi/2 - polarAngle, (N, 1))
    R = np.tile(r.reshape(-1, 1), (1, nshot*nseg))

    x, y, z = sph2cart(azimuthal, polar, R)

    x = x.reshape(N, nseg, nshot)
    y = y.reshape(N, nseg, nshot)
    z = z.reshape(N, nseg, nshot)

    if flagPlot:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for shot in range(20):  # 201:nshot
            for seg in range(nseg):
                ax.plot3D(x[:, seg, shot], y[:, seg, shot], z[:, seg, shot])
                ax.set_title('N = {}  nseg = {}  nshot = {}'.format(N, nseg, nshot))
                ax.set_xlabel('X')
                ax.set_ylabel('Y')
                ax.set_zlabel('Z')
                ax.set_xlim(-0.5, 0.5)
                ax.set_ylim(-0.5, 0.5)
                ax.set_zlim(-0.5, 0.5)
                plt.pause(0.1)
        plt.show()
    
    return x, y, z, polarAngle, azimuthalAngle

def get_np_array_SiemensrawData(rawDataPath):
    multi_twix: list = read_twix(rawDataPath, include_scans=[-1], keep_syncdata_and_acqend=False)
    datas = [mdb.data for mdb in multi_twix[0]['mdb']]
    return np.array(datas)

def get_np_array_SiemensTimestamps(rawDataPath):
    multi_twix: list = read_twix(rawDataPath, include_scans=[-1], keep_syncdata_and_acqend=False)
    timestamps = [mdb.mdh.TimeStamp for mdb in multi_twix[0]['mdb']]
    pmuTimestamps = [mdb.mdh.PMUTimeStamp for mdb in multi_twix[0]['mdb']]
    return np.array(timestamps), np.array(pmuTimestamps)

def get_np_array_SiemensPmuTimestamps(rawDataPath):
    multi_twix: list = read_twix(rawDataPath, include_scans=[-1], keep_syncdata_and_acqend=False)
    pmuTimestamps = [mdb.mdh.PMUTimeStamp for mdb in multi_twix[0]['mdb']]
    return np.array(pmuTimestamps)

def loadNpArrayFormatted(path):
    """ PLEASE CHECK CORRECTNESS OF THIS FUNCTION if you want to use it, it's just a draft to show how to load the data can be combined with the phyllotaxis function
    """
    arr = get_np_array_SiemensrawData(path)
    
    # Assert that shape of reshaped_array first dimention is a mupltiple of 22
    assert arr.shape[0] % 22 == 0
    N = 480
    nseg = 22
    nshot = arr.shape[0]//22
    flagSelfNav = True
    flagPlot = False
    # Generate the phyllotaxis trajectory with 22, rawData.shape[0]//22
    x, y, z, _ , _ =  computePhyllotaxis(N,nseg, nshot, flagSelfNav, flagPlot)

    kSpaceCoordinates = np.stack((x,y,z), axis=-1)
    kSpaceCoordinatesFlat = kSpaceCoordinates.transpose(0, 2, 1, 3).reshape(480, -1, 3)
    kSpaceCoordinatesFlat.shape
    mergedArray = np.concatenate((kSpaceCoordinatesFlat, arr.transpose(2, 0, 1)), axis=-1).transpose(1, 2, 0)
    
    pmu = get_np_array_SiemensPmuTimestamps(path)
    pmu = pmu[:, np.newaxis, np.newaxis]
    # Repeat timestamps along the second axis to match the shape of datas
    pmu = np.repeat(pmu, mergedArray.shape[1], axis=1)
    # Append timestamps to the datas array along the last axis
    arr_with_pmu = np.concatenate((mergedArray, pmu), axis=-1)
    arr_with_pmu = arr_with_pmu.reshape(arr_with_pmu.shape[0], -1)
    
    return arr_with_pmu