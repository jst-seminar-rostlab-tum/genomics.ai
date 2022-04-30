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
    description: 'Test',
    profilePictureURL: 'https://www.hzdr.de/db/Pic?pOid=55058',
    backgroundPictureURL: 'https://upload.wikimedia.org/wikipedia/commons/thumb/2/22/Helmholtz_Zentrum_M%C3%BCnchen.jpg/2560px-Helmholtz_Zentrum_M%C3%BCnchen.jpg',
    adminIds: [1,4,5],
    memberIds: [2, 3],
  },
  {
    id: 2,
    name: 'Technische Universität München',
    country: 'Germany',
    profilePictureURL: 'https://upload.wikimedia.org/wikipedia/commons/b/ba/Tum_logo.gif',
    backgroundPictureURL: 'https://www.in.tum.de/fileadmin/_processed_/5/5/csm_2006_1015Bild0136_9dc504e910.jpg',
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

export async function removeMemberFromInstitution(institutionId, memberId) {
  const institution = testInstitutions.find((i) => i.id === institutionId);
  if (!institution) {
    throw new Error('Institution not found');
  }
  institution.memberIds = institution.memberIds.filter((id) => id !== memberId);
}