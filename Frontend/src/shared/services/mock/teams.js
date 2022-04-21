const mockLeftIds = [];
let runningId = 1;

export async function leaveTeam(team) {
  mockLeftIds.push(team.id);
}

export async function createTeam(name, description) {
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
  };
}

export async function getTeam(id) {
  if (parseInt(id, 10) === 1) {
    return {
      id: 1,
      name: 'Biotechnology Team',
      description: 'Biotechnology Team bla bla',
      adminIds: [1],
      invitedMemberIds: [],
      memberIds: [1, 2],
      visibility: 'public',
      institutionId: 2,
    };
  }
  throw Error('TeamNotFound (mock)');
}

export default async function queryMyTeams() {
  return [
    {
      id: 1,
      name: 'Biotechnology Team',
      description: 'Biotechnology Team bla bla',
      adminIds: [1],
      invitedMemberIds: [],
      memberIds: [1, 2],
      visibility: 'public',
      institutionId: 2,
    },
  ].filter((team) => !mockLeftIds.includes(team.id));
}
