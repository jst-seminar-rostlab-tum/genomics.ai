const mockLeftIds = [];

export async function leaveTeam(team) {
  mockLeftIds.push(team.id);
}

export async function getTeam(id) {
  if (parseInt(id, 10) === 1) {
    return {
      id: 1,
      name: 'Biotechnology Team',
      description: 'Biotechnology Team bla bla',
      adminIds: [1],
      invitedMemberIds: [],
      memberIds: [2],
      visibility: 'public',
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
      memberIds: [2],
      visibility: 'public',
    },
  ].filter((team) => !mockLeftIds.includes(team.id));
}
