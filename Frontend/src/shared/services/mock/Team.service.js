import MemberService from './Member.service';
import ProfileService from './Profile.service';

let mockTeams = [
  {
    id: '1',
    name: 'Biotechnology Team',
    description: 'Biotechnology Team bla bla',
    adminIds: [6, 7],
    invitedMemberIds: [],
    memberIds: [1, 2, 3, 4, 5],
    projects: [],
    visibility: 'public',
    institutionId: 1,
  },
];
let runningId = 1;

const TeamService = {
  async leaveTeam(team) {
    mockTeams = mockTeams.filter((t) => t.id !== team.id);
  },

  async createTeam(name, description, institutionId) {
    // fake effect
    await new Promise((resolve) => setTimeout(resolve, 1000));
    runningId += 1;
    mockTeams.push(
      {
        id: runningId.toString(),
        name,
        country: null,
        description,
        avatarUrl: null,
        backgroundPictureURL: null,
        adminIds: [1], // TODO: make sure that the backend puts my user ID here
        memberIds: [],
        visibility: 'private',
        institutionId,
      },
    );
    return mockTeams[mockTeams.length - 1];
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
    return preparedTeams.sort((a, b) => (`${a.title}`).localeCompare(b.title));
  },

};

export default TeamService;
