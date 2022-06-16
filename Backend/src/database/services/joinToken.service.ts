import { AddInstitutionJoinTokenDTO, AddTeamJoinTokenDTO } from "../dtos/joinToken.dto";
import { joinTokenModel, JoinTarget, IJoinToken } from "../models/join_token";

export default class JoinTokenService {
  static async addInstitutionJoinToken(token: AddInstitutionJoinTokenDTO): Promise<IJoinToken> {
    return await joinTokenModel.create({
      userId: token.userId,
      institutionId: token.institutionId,
      target: JoinTarget.INSTITUTION,
    });
  }

  static async addTeamJoinToken(token: AddTeamJoinTokenDTO): Promise<IJoinToken> {
    return await joinTokenModel.create({
      userId: token.userId,
      teamId: token.teamId,
      target: JoinTarget.TEAM,
    });
  }

  static async getTokenData(token: string): Promise<IJoinToken> {
    return await joinTokenModel.findOne({ token });
  }

  static async removeToken(token: string): Promise<any> {
    return await joinTokenModel.deleteOne({ token });
  }
}
