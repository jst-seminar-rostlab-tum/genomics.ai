const mockLeftIds = [];

export async function leaveProject(project) {
  mockLeftIds.push(project.id);
}

export default async function queryMyProjects() {
  return [
    {
      id: 1,
      name: 'Biotechnology Project',
      description: 'Biotechnology Project bla bla',
      adminIds: [1],
      invitedMemberIds: [],
      memberIds: [2],
      visibility: 'public',
    },
  ].filter((project) => !mockLeftIds.includes(project.id));
}
