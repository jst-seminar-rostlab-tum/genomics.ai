const mockProjects = [
  {
    "_id": "626f20715c2633c6573652d3",
    "owner": "626bdb1ed76c8b968a50f833",
    "name": "myproject",
    "modelId": "626ea2361d7d1a27de465b5f",
    "atlasId": "626ea3311d7d1a27de465b63",
    "fileName": "pancreas_query.h5ad",
    "fileSize": 21892149,
    "uploadDate": "2022-05-02T00:06:09.219Z",
    "status": "DONE",
    "resultSize": -1,
    "__v": 0,
    "uploadId": "ABPnzm64TXY8lmYTzYiCOoX4J5S-VsKxK26LilPb3aex0I6LOtNXBVyU8XnT88ns76ZBY3cz",
    "teamId": "626edd8b75993ffbfb5ae41d",
  },
  {
    "_id": "626fc5d175993ffbfb5ae46c",
    "owner": "626bdb1ed76c8b968a50f833",
    "name": "This is Project",
    "modelId": "626ea2361d7d1a27de465b5f",
    "atlasId": "626ea3311d7d1a27de465b63",
    "fileName": "pbmc_query.h5ad",
    "fileSize": 21892149,
    "uploadDate": "2022-05-02T00:06:09.219Z",
    "status": "PROCESSING_PENDING",
    "resultSize": -1,
    "__v": 0,
    "uploadId": "ABPnzm64TXY8lmYTzYiCOoX4J5S-VsKxK26LilPb3aex0I6LOtNXBVyU8XnT88ns76ZBY3cz",
    "teamId": "626edd8b75993ffbfb5ae41d",
  },
];

const ProjectService = {
  async getProjects() {
    return mockProjects;
  },

  async getTeamProjects(teamId) {
    return mockProjects;
  },

  getOwnProjects: async () => ProjectService.getProjects(),

  async getProject(id) {
    const project = mockProjects.find((p) => String(p._id) === String(id));
    return project ?? {
      _id: id,
      owner: 'mock',
      uploadId: '12345-abcde',
      atlasId: '626ea3311d7d1a27de465b63',
      fileName: 'test1.h5ad',
      location: `./testData/test_file${id}.csv`,
      fileSize: 10000000,
      modelId: '626ea2361d7d1a27de465b5e',
      name: `Test project ${id}`,
      resultSize: '10000',
      status: 'DONE',
      uploadDate: '2022-04-30T106:00:00.000Z',
    };
  },
};

export default ProjectService;
