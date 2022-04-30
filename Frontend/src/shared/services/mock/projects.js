const projects = [
  {
    _id: 1,
    owner: 'mock',
    uploadId: '12345-abcde',
    atlasId: 3,
    fileName: 'test1.h5ad',
    location: './testData/test_file1.csv',
    fileSize: 10000000,
    modelId: 2,
    name: 'test project 01',
    resultSize: '10000',
    status: 'DONE',
    uploadDate: '2022-04-30T106:00:00.000Z',
  },
  {
    _id: 2,
    owner: 'mock',
    uploadId: 'abcde-12345',
    atlasId: 1,
    fileName: 'test2.h5ad',
    location: './testData/test_file2.csv',
    fileSize: 10000000,
    modelId: 4,
    name: 'test project 02',
    resultSize: '10000',
    status: 'DONE',
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
    _id: 1,
    name: 'Atlas 1',
    previewPictureURL: './previewImages/atlas.jpg',
    modalities: ['RNA', 'ADT'],
    numberOfCells: 5348812,
    species: ['Human'],
    compatibleModels: [1, 2, 3],
  },
  {
    _id: 2,
    name: 'Atlas 2',
    previewPictureURL: './previewImages/atlas.jpg',
    modalities: ['RNA', 'ADT'],
    numberOfCells: 15981239,
    species: ['Mouse'],
    compatibleModels: [1, 2, 3],
  },
  {
    _id: 3,
    name: 'Atlas 3',
    previewPictureURL: './previewImages/atlas.jpg',
    modalities: ['RNA'],
    numberOfCells: 693493223,
    species: ['Human'],
    compatibleModels: [1, 2],
  },
  {
    _id: 4,
    name: 'Atlas 4',
    previewPictureURL: './previewImages/atlas.jpg',
    modalities: ['RNA'],
    numberOfCells: 13453230,
    species: ['Human'],
    compatibleModels: [1, 2, 3],
  },
  {
    _id: 5,
    name: 'Atlas 5',
    previewPictureURL: './previewImages/atlas.jpg',
    modalities: ['RNA'],
    numberOfCells: 86534143,
    species: ['Human'],
    compatibleModels: [1, 3],
  },
];

const models = [
  {
    _id: 1,
    name: 'Model 1',
    description: 'Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Eget nunc lobortis mattis aliquam faucibus purus in. Interdum velit laoreet id donec ultrices tincidunt arcu. Proin sagittis nisl rhoncus mattis rhoncus urna neque viverra. Interdum posuere lorem ipsum dolor. Tempus iaculis urna id volutpat lacus laoreet non curabitur. Suscipit tellus mauris a diam maecenas sed enim ut. Elementum nibh tellus molestie nunc non. Malesuada fames ac turpis egestas sed.',
    requirements: ['h5ad file', 'labeled data'],
  },
  {
    _id: 2,
    name: 'Model 2',
    description: 'Lorem ipsum dolor sit amet.',
    requirements: ['h5', 'labeled data'],
  },
  {
    _id: 3,
    name: 'Model 3',
    description: 'Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Eget nunc lobortis mattis aliquam faucibus purus in. Interdum velit laoreet id donec ultrices tincidunt arcu. Proin sagittis nisl rhoncus mattis rhoncus urna neque viverra.',
    requirements: ['h5', 'unlabeled data'],
  },
];

export default async function getProjects() {
  return projects;
}

export async function getProject(id) {
  await new Promise((resolve) => setTimeout(resolve, 1000));
  return projects.find((project) => project._id === id);
}

export async function addProjectToTeam(teamId, projectId) {
  throw new Error('400');
}

export async function getOwnTeams() {
  return ownTeams;
}

export async function getAtlas(id) {
  return atlases.find((atlas) => atlas._id === id);
}

export async function getModels() {
  return models;
}

export async function getAtlases() {
  return atlases;
}