import { ObjectId } from "mongoose";


export interface AddInstitutionJoinTokenDTO {
  userId: ObjectId;
  institutionId: ObjectId;
}
export interface AddTeamJoinTokenDTO {
  userId: ObjectId;
  teamId: ObjectId;
}
