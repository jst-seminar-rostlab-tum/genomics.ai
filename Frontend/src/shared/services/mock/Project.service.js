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
};

export default ProjectService;
