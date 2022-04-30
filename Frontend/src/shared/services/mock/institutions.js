const mockLeftIds = [];
let runningId = 3;

export async function createInstitution(name, description) {
  // fake effect
  await new Promise((resolve) => setTimeout(resolve, 1000));
  runningId += 1;
  return {
    id: runningId,
    name,
    country: null,
    description,
    profilePictureURL: null,
    backgroundPictureURL: null,
    adminIds: [1], // TODO: make sure that the backend puts my user ID here
    memberIds: [],
  };
}

const testInstitutions = [
  {
    id: 1,
    name: 'Helmholtz Institute',
    country: 'Germany',
    profilePictureURL: 'https://www.hzdr.de/db/Pic?pOid=55058',
    backgroundPictureURL: 'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAAAXNSR0IArs4c6QAAAA1JREFUGFdjsC9c9x8ABK8CXkzrf1cAAAAASUVORK5CYII=',
    adminIds: [],
    memberIds: [],
  },
  {
    id: 2,
    name: 'Technische Universität München',
    country: 'Germany',
    profilePictureURL: 'https://scalings.eu/wp-content/uploads/2019/11/tum-logo.png',
    backgroundPictureURL: 'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAAAXNSR0IArs4c6QAAAA1JREFUGFdjsC9c9x8ABK8CXkzrf1cAAAAASUVORK5CYII=',
    adminIds: [1],
    memberIds: [2, 3, 4, 5],
  },
  {
    id: 3,
    name: 'Rostlab',
    country: 'Germany',
    profilePictureURL: 'https://avatars.githubusercontent.com/u/4093405?s=200&v=4',
    backgroundPictureURL: 'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAAAXNSR0IArs4c6QAAAA1JREFUGFdjsC9c9x8ABK8CXkzrf1cAAAAASUVORK5CYII=',
    adminIds: [1],
    memberIds: [2, 3, 4, 5],
  },
];

export async function leaveInstitution(institution) {
  mockLeftIds.push(institution.id);
}

export default async function queryMyInstitutions() {
  return testInstitutions.filter((institution) => !mockLeftIds.includes(institution.id));
}

export async function queryIsAdminInstitutions(userId) {
  return testInstitutions.filter((institution) => institution.adminIds.includes(userId));
}

export async function getInstitution(id) {
  return testInstitutions.find((institution) => institution.id === id);
}
