import MemberService from '../Member.service';
import ProfileService from '../Profile.service';

let mockTeams = [
  {
    id: '1',
    name: 'Biotechnology Team',
    description: 'The Biotechnology Team works in integration of natural sciences and engineering sciences in order to achieve the application of organisms, cells, parts thereof and molecular analogues for products and services.',
    adminIds: [6, 7],
    invitedMemberIds: [],
    memberIds: ['626bdb1ed76c8b968a50f833', 2, 3, 4, 5],
    projects: [],
    visibility: 'public',
    institutionId: '1',
  },
  {
    id: '2',
    name: 'Genomics Team',
    description: 'The Genomics Team works in field of biology focusing on the structure, function, evolution, mapping, and editing of genomes.',
    adminIds: ['626bdb1ed76c8b968a50f833'],
    invitedMemberIds: [],
    memberIds: [2, 3, 4, 5],
    projects: [],
    visibility: 'public',
    institutionId: '3',
  },
  {
    id: '3',
    name: 'Heinig Lab',
    description: 'Heinig Lab bla bla',
    adminIds: ['626bdb1ed76c8b968a50f833'],
    invitedMemberIds: [],
    memberIds: [3, 4, 5],
    projects: [],
    visibility: 'private',
    institutionId: '3',
  },
];
let runningId = 2;

const TeamService = {
  async leaveTeam(team) {
    mockTeams = mockTeams.filter((t) => t.id !== team.id);
  },

  async createTeam(name, description, institutionId) {
    // fake effect
    await new Promise((resolve) => setTimeout(resolve, 1000));
    runningId += 1;
    const user = await ProfileService.getProfile();
    const newTeam = {
      id: runningId.toString(),
      name,
      country: null,
      description,
      backgroundPictureURL: null,
      adminIds: [user.id], // TODO: make sure that the backend puts my user ID here
      memberIds: [],
      invitedMemberIds: [],
      visibility: 'private',
      institutionId,
    };
    mockTeams = [...mockTeams, newTeam];
    return newTeam;
  },

  async removeMemberFromTeam(teamId, memberId) {
    const team = mockTeams.find((t) => t.id === teamId);
    if (!team) {
      throw new Error('Team not found');
    }
    team.memberIds = team.memberIds.filter((id) => id !== memberId);
  },

  async getTeam(id) {
    const team = mockTeams.find((t) => t.id === id);
    if (team == null) {
      throw new Error(`Could not find team with id ${id}`);
    }
    return team;
  },

  async getInstitutionTeams(institutionId) {
    return mockTeams.filter((team) => team.institutionId === institutionId);
  },

  async getMembers(teamId) {
    const team = await TeamService.getTeam(teamId);
    return Promise.all([
      ...team.adminIds,
      ...team.memberIds,
    ].map(MemberService.getMember));
  },

  async getMyTeams() {
    return mockTeams;
  },

  async getMyAdminTeams() {
    const myTeams = await TeamService.getMyTeams();
    const user = await ProfileService.getProfile();
    return myTeams.filter((t) => t.adminIds.includes(user.id));
  },
  getTeams: async (params) => {
    let preparedTeams = mockTeams.map(
      (team) => ({ ...team, title: team.name }
      ),
    );
    if (params.keyword) {
      preparedTeams = preparedTeams.filter(
        (team) => team.name.toLowerCase().includes(params.keyword.toLowerCase()),
      );
    }
    if (params.visibility) {
      preparedTeams = preparedTeams.filter(
        (team) => team.visibility.toUpperCase().includes(params.visibility),
      );
    }

    if (!params.sortBy || params.sortBy === 'name') {
      preparedTeams = preparedTeams.sort((a, b) => (`${a.name}`).localeCompare(b.name));
    }
    return preparedTeams;
  },

};

export default TeamService;
