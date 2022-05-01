import React, { useState, useEffect } from 'react';
import HeaderView from 'components/general/HeaderView';
import InstitutionList from 'components/institutions/InstitutionList';
import TabPanel from 'components/TabPanel'
import { useAuth } from 'shared/context/authContext';
import Box from '@mui/material/Box';
import Tabs from '@mui/material/Tabs';
import Tab from '@mui/material/Tab';
import Grid from '@mui/material/Grid';
import ProfileImage from 'components/ProfileImage';
import TeamCard from 'components/teams/TeamCard';

const testInstitutions = [
    {
        id: 5,
        name: 'Helmholtz Institute',
        country: 'Germany',
        profilePictureURL: 'https://www.hzdr.de/db/Pic?pOid=55058',
        backgroundPictureURL: 'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAAAXNSR0IArs4c6QAAAA1JREFUGFdjsC9c9x8ABK8CXkzrf1cAAAAASUVORK5CYII=',
        adminIds: [],
        memberIds: [],
    },
    {
        id: 6,
        name: 'Technische Universität München',
        country: 'Germany',
        profilePictureURL: 'https://scalings.eu/wp-content/uploads/2019/11/tum-logo.png',
        backgroundPictureURL: 'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAAAXNSR0IArs4c6QAAAA1JREFUGFdjsC9c9x8ABK8CXkzrf1cAAAAASUVORK5CYII=',
        adminIds: [1],
        memberIds: [],
    },
    {
        id: 7,
        name: 'Rostlab',
        country: 'Germany',
        profilePictureURL: 'https://avatars.githubusercontent.com/u/4093405?s=200&v=4',
        backgroundPictureURL: 'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAAAXNSR0IArs4c6QAAAA1JREFUGFdjsC9c9x8ABK8CXkzrf1cAAAAASUVORK5CYII=',
        adminIds: [1],
        memberIds: [],
    },
];

let mockTeams = [
    {
        id: 1,
        name: 'Biotechnology Team',
        description: 'Biotechnology Team bla bla',
        adminIds: [6, 7],
        invitedMemberIds: [],
        memberIds: [1, 2, 3, 4, 5],
        visibility: 'public',
        institutionId: 2,
    },
];

function UserProfile() {
    const [user] = useAuth()

    const [value, setValue] = React.useState(0);

    const handleChange = (event, newValue) => {
        setValue(newValue);
    };

    const [institutions, setInstitutions] = useState(testInstitutions);
    const [teams, setTeams] = useState(mockTeams);

    return (
        <>
            <HeaderView
                title="Profile"
            >
                <Grid container justifyContent="center" alignItems="center" spacing={3}>
                    <Grid item>
                        <ProfileImage sizePixels={160} />
                    </Grid>
                    <Grid item>
                        <h2>{user.firstName} {user.lastName}</h2>
                        <span>{user.email}</span>
                    </Grid>
                </Grid>
                <Box sx={{ width: '100%' }}>
                    <Box sx={{ borderBottom: 1, borderColor: 'divider' }}>
                        <Tabs value={value} onChange={handleChange} aria-label="basic tabs example">
                            <Tab label="Institutions" {...a11yProps(0)} />
                            <Tab label="Teams" {...a11yProps(1)} />
                        </Tabs>
                    </Box>
                    <TabPanel value={value} index={0}>
                        <InstitutionList
                            isLoading={false}
                            institutions={institutions}
                        />
                    </TabPanel>
                    <TabPanel value={value} index={1}>
                        {teams.length === 0 ? 'No teams.' : ''}
                        {teams.map((team) => (
                            <div key={team.id}>
                                <TeamCard
                                    team={team}
                                    onLeft={(t) => console.log("")}
                                />
                            </div>
                        ))}

                    </TabPanel>
                </Box>
            </HeaderView>
        </>
    );
}

function a11yProps(index) {
    return {
        id: `simple-tab-${index}`,
        'aria-controls': `simple-tabpanel-${index}`,
    };
}

export default UserProfile;
