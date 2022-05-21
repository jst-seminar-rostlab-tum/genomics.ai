import React, { useState, useEffect } from 'react';
import { useParams } from 'react-router-dom';
import HeaderView from 'components/general/HeaderView';
import TeamProjectList from 'components/teams/detail/TeamProjectList';
import TeamMemberList from 'components/teams/detail/TeamMemberList';
import TeamAdminHeaderRight from 'components/teams/detail/TeamAdminHeaderRight';
import TeamUserHeaderRight from 'components/teams/detail/TeamUserHeaderRight';
import TeamHeaderOptions from 'components/teams/detail/TeamHeaderOptions';
import TeamService from 'shared/services/Team.service';
import InstitutionService from 'shared/services/Institution.service';
import TeamInviteButton from 'components/teams/detail/TeamInviteButton';
import TextField from '@mui/material/TextField';
import { useAuth } from 'shared/context/authContext';

export default function TeamPage() {
  const { id } = useParams();
  const [team, setTeam] = useState({});
  const [user] = useAuth();
  const [institution, setInstitution] = useState({});

  const handleDescriptionChange = (event) => {
    setTeam({
      ...team,
      description: event.target.value,
    });
  };

  useEffect(() => {
    if (id == null) return;
    TeamService.getTeam(id)
      .then(setTeam)
      .catch((ignored) => { console.error(ignored); });
  }, [setTeam]);

  // Institution may be undefined
  useEffect(() => {
    if (team.institutionId) {
      InstitutionService.getInstitution(team.institutionId)
        .then((newInstitution) => setInstitution(newInstitution));
    }
  }, [team]);

  const isAdmin = team.adminIds ? team.adminIds.includes(user._id) : false;

  if (!team) return <></>;

  return (
    <HeaderView
      title={team.name}
      rightOfTitle={(
        <TeamHeaderOptions
          team={team}
          isAdmin={isAdmin}
          institution={institution}
        />
      )}
      replaceHeaderRight={
        (isAdmin && <TeamAdminHeaderRight team={team} setTeam={setTeam} />)
        || <TeamUserHeaderRight institution={institution} team={team} user={user} />
      }
    >
      <br />
      <section>
        <h2>Description</h2>
        <hr />
        <TextField
          id="description"
          multiline
          minRows={3}
          maxRows={5}
          value={team.description}
          InputProps={{
            readOnly: !isAdmin,
          }}
          fullWidth
          onChange={handleDescriptionChange}
          variant="standard"
        />
      </section>
      <section>
        <h2>GeneMapper</h2>
        <hr />
        <TeamProjectList team={team} forPart="geneMapper" />
      </section>
      <section>
        <h2>GeneCruncher</h2>
        <hr />
        <TeamProjectList team={team} forPart="geneCruncher" />
      </section>
      <section>
        <h2>Members</h2>
        <hr />
        <TeamMemberList
          team={team}
          onMemberRemoved={(_team, member) => {
            setTeam({
              ...team,
              memberIds: team.memberIds.filter((mId) => mId !== member.id),
              adminIds: team.adminIds.filter((aId) => aId !== member.id),
            });
          }}
          onMakeAdmin={(_team, member) => {
            setTeam({
              ...team,
              memberIds: team.memberIds.filter((mId) => mId !== member.id),
              adminIds: [...team.adminIds, member.id],
            });
          }}
          onRemoveAdmin={(_team, member) => {
            setTeam({
              ...team,
              memberIds: [...team.memberIds, member.id],
              adminIds: team.adminIds.filter((aId) => aId !== member.id),
            });
          }}
        />
      </section>
      {isAdmin && <TeamInviteButton team={team} />}
    </HeaderView>
  );
}
