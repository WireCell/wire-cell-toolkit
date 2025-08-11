import torch, numpy as np, matplotlib.pyplot as plt
from argparse import ArgumentParser as ap

if __name__ == '__main__':

    parser = ap()
    parser.add_argument('--no-activities', help="Don't show activities", action='store_true')
    parser.add_argument('--no-blobs', help="Don't show blobs", action='store_true')
    parser.add_argument('--max-blobs', type=int, help='Max number of blobs to draw', default=10)
    parser.add_argument('--skip-views', type=int, nargs='*')
    args = parser.parse_args()

    points = torch.load('random_points.zip').numpy()
    blobs = torch.load('solved_blobs.zip').numpy()
    views = torch.load('views.zip').numpy()
    print('Views')
    print(views)
    #Scatter the points
    n=1000

    dirs = views[:, 1] - views[:, 0]
    print('Dirs')
    print(dirs)

    # Define the 90-degree counter-clockwise rotation matrix
    rotation_matrix_ccw = np.array(
        [[0, -1],
        [1,  0]]
    )
    rot = np.zeros_like(dirs)
    for i in range(dirs.shape[0]): rot[i] = np.dot(rotation_matrix_ccw, dirs[i])

    print(rot)

    colors = ['k', 'grey', 'r', 'c', 'm']

    if not args.no_activities:
        for iview in range(2, 5):
            if args.skip_views is not None and iview in args.skip_views: continue
            active = torch.load(f'activities{iview}.zip').nonzero().numpy()
            active_locs = views[iview, 0] + dirs[iview]*active
            # print(active_locs)
            for iactive in range(active_locs.shape[0]):
                plt.axline(active_locs[iactive], active_locs[iactive]+rot[iview], c=colors[iview], linestyle='--')

    used_blobs = [[] for i in range(views.shape[0])]


    if not args.no_blobs:
        #iterate over the blobs
        print(f'Found {blobs.shape[0]} blobs')
        for i in range(min([args.max_blobs, blobs.shape[0]])):
            blob = blobs[i]
            # print('Blob\n', blob)
            blob_starts = views[:,0] + dirs*blob[:,0].reshape(-1, 1)
            blob_ends = views[:,0] + dirs*blob[:,1].reshape(-1, 1)
            # print('Blob starts\n', blob_starts)
            # print('Blob ends\n', blob_ends)

            for j in range(blob.shape[0]):
                if args.skip_views is not None and j in args.skip_views: continue
                if blob[j].tolist() in used_blobs[j]: continue
                plt.axline(blob_starts[j], blob_starts[j]+rot[j], c=colors[j])
                plt.axline(blob_ends[j], blob_ends[j]+rot[j], c=colors[j])
                used_blobs[j].append(blob[j].tolist())

    plt.scatter(points[:n, :,0], points[:n, :,1], s=10, c='k')
    plt.xlim(-50, 4050)
    plt.ylim(-50, 4050)
    plt.show()