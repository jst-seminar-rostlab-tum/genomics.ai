const projects = [
  {
    _id: 1,
    owner: 'mock',
    uploadId: '12345-abcde',
    atlasId: '626ea3311d7d1a27de465b63',
    fileName: 'test1.h5ad',
    location: './testData/test_file1.csv',
    fileSize: 10000000,
    modelId: '626ea2361d7d1a27de465b5e',
    name: 'Demo project 01',
    resultSize: '10000',
    status: 'DONE',
    uploadDate: '2022-04-30T106:00:00.000Z',
  },
  {
    _id: 2,
    owner: 'mock',
    uploadId: 'abcde-12345',
    atlasId: '626ea3311d7d1a27de465b64',
    fileName: 'test2.h5ad',
    location: './testData/test_file2.csv',
    fileSize: 10000000,
    modelId: '626ea2361d7d1a27de465b60',
    name: 'Another sample project',
    resultSize: '10000',
    status: 'DONE',
    uploadDate: '2022-04-30T106:00:00.000Z',
  },
  {
    _id: 3,
    owner: 'mock',
    uploadId: 'abcde-12345',
    atlasId: '626ea3311d7d1a27de465b64',
    fileName: 'test3.h5ad',
    location: './testData/test_file3.csv',
    fileSize: 10000000,
    modelId: '626ea2361d7d1a27de465b60',
    name: 'Test project 03',
    resultSize: '10000',
    status: 'PROCESSING_PENDING',
    uploadDate: '2022-04-30T106:00:00.000Z',
  },
  {
    _id: 4,
    owner: 'mock',
    uploadId: 'abcde-12345',
    atlasId: '626ea3311d7d1a27de465b63',
    fileName: 'test3.h5ad',
    location: './testData/test_file4.csv',
    fileSize: 10000000,
    modelId: '626ea2361d7d1a27de465b5f',
    name: 'Sample project',
    resultSize: '10000',
    status: 'ABORTED',
    uploadDate: '2022-04-30T106:00:00.000Z',
  },
];

const ownTeams = [
  {
    _id: 1,
    title: 'test team 01',
  },
  {
    _id: 2,
    title: 'test team 02',
  },
];

const atlases = [
  {
    _id: '626ea3311d7d1a27de465b63',
    name: 'Pancreas',
    previewPictureURL: 'https://images-na.ssl-images-amazon.com/images/I/81+6yFP1sQL.jpg',
    modalities: ['No'],
    numberOfCells: 10293,
    species: ['Human'],
    compatibleModels: ['scVI', 'scanVI'],
  },
  {
    _id: '626ea3311d7d1a27de465b64',
    name: 'PBMC',
    previewPictureURL: 'https://images-na.ssl-images-amazon.com/images/I/81+6yFP1sQL.jpg',
    modalities: ['CITE-seq: CD3_TotalSeqB, CD4_TotalSeqB, CD8a_TotalSeqB, CD14_TotalSeqB, CD15_TotalSeqB, CD16_TotalSeqB, CD56_TotalSeqB, CD19_TotalSeqB, CD25_TotalSeqB, CD45RA_TotalSeqB, CD45RO_TotalSeqB, PD-1_TotalSeqB, TIGIT_TotalSeqB, CD127_TotalSeqB'],
    numberOfCells: 10849,
    species: ['Human'],
    compatibleModels: ['totalVI'],
  }];

const models = [
  {
    _id: '626ea2361d7d1a27de465b5e',
    name: 'scANVI',
    description: 'The scANVI model lorem impsum blablabla lorem impsum blablabla lorem impsum blablabla lorem impsum blablabla lorem impsum blablabla lorem impsum blablabla',
  },
  {
    _id: '626ea2361d7d1a27de465b5f',
    name: 'scVI',
    description: 'The scVI model lorem impsum blablabla lorem impsum blablabla lorem impsum blablabla lorem impsum blablabla lorem impsum blablabla lorem impsum blablabla',
  },
  {
    _id: '626ea2361d7d1a27de465b60',
    name: 'totalVI',
    description: 'The totalVI model lorem impsum blablabla lorem impsum blablabla lorem impsum blablabla lorem impsum blablabla lorem impsum blablabla lorem impsum blablabla',
  },
];

const datasets = [
  {
    _id: 1,
    name: 'Hao and Hao et al, bioRvix 2020',
    status: 'DONE',
  },
  {
    _id: 2,
    name: 'Keller, Lennart et al, FastGenomics.ai 2022',
    status: 'IN PROGRESS',
  },
  {
    _id: 3,
    name: 'Hübsch und Hübner et al, ToTheMoon.crypto 2021',
    status: 'UPLOAD FAILED',
  },
  {
    _id: 4,
    name: 'Xu Fan Lu et al, BINGO 2022',
    status: 'DONE',
  },
  {
    _id: 5,
    name: 'Christopher, Evan et al, TUM.UNICORN 2023',
    status: 'IN PROGRESS',
  },
];

const ProjectMock = {
  getProjects: async () => projects,
  getProject: async (id) => projects.find((project) => project._id === Number(id)),

  addProjectToTeam: async (teamId, projectId) => { throw new Error('400'); },
  getOwnTeams: async () => ownTeams,

  getAtlases: async () => atlases,
  getAtlas: async (id) => atlases.find((atlas) => atlas._id === Number(id)),

  getModel: async (id) => models.find((model) => model._id === Number(id)),
  getModels: async () => models,

  getDataset: async (id) => models.find((dataset) => dataset._id === Number(id)),
  getDatasets: async () => datasets,

  deleteProject: async (id) => {
    const i = projects.findIndex((project) => project._id === Number(id));
    return i > -1 ? projects.splice(i, 1) : false;
  },
};

export default ProjectMock;
